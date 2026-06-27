// ── src/fastEditEdges.cpp ───────────────────────────────────────────────────
#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Lift floating-point products back over an integer they should have reached:
// thresh * norm_len can land at e.g. 2.9999999996 when the true value is 3, and
// a bare (int) cast would drop a pair whose normalized distance exactly equals
// the threshold (the inclusive `fd <= thresh` boundary). The epsilon dwarfs the
// FP rounding error (~norm_len * 2e-16) while never reaching the next integer.
static const double THRESH_EPS = 1e-9;

// ============================================================================
//  Deletion-neighborhood (SymSpell) helpers
// ============================================================================
// All strings reachable by deleting up to `k` characters. For Levenshtein,
// edit(a,b) <= k implies their deletion neighborhoods intersect, so colliding
// reps are exact candidates (no false negatives) verified later by banded DP.
// A variant of length L-d is always reached with exactly d deletions, so the
// "insert-then-recurse only on first sight" prune is loss-free.
static void gen_deletes_rec(const std::string& cur, int rem,
                            std::unordered_set<std::string>& out) {
  if (rem == 0) return;
  for (size_t i = 0; i < cur.size(); ++i) {
    std::string nxt;
    nxt.reserve(cur.size() - 1);
    nxt.append(cur, 0, i);
    nxt.append(cur, i + 1, std::string::npos);
    if (out.insert(nxt).second) gen_deletes_rec(nxt, rem - 1, out);
  }
}
static void gen_deletes(const std::string& s, int k,
                        std::unordered_set<std::string>& out) {
  out.clear();
  out.insert(s);
  if (k > 0) gen_deletes_rec(s, k, out);
}
// Rough upper bound on variants per sequence: sum_{j=0}^{k} C(L, j), capped.
static long deletion_variant_bound(int L, int k) {
  long tot = 1, c = 1;
  for (int j = 1; j <= k; ++j) {
    c = c * (L - j + 1) / j;
    tot += c;
    if (tot > 1000000L) return tot;
  }
  return tot;
}

// ============================================================================
//  Data Structures & Buffers (For Speed)
// ============================================================================

struct SubstMatrix {
  std::vector<std::vector<int>> scores;
  std::vector<int> char_map; // Faster than unordered_map for lookups
  int gap_open;
  int gap_extend;
  
  SubstMatrix() : char_map(256, -1) {} // Map ASCII chars to indices
  
  inline int get_score(unsigned char a, unsigned char b) const {
    int i = char_map[a];
    int j = char_map[b];
    if (i < 0 || j < 0) return -4; 
    return scores[i][j];
  }
};

// Thread-local workspace to prevent re-allocation in loops
struct DPWorkspace {
  std::vector<int> v1;
  std::vector<std::vector<int>> mat;
  std::vector<int> last_row_tracker; // For Damerau
  
  void resize_mat(int n, int m) {
    if ((int)mat.size() <= n) mat.resize(n + 5);
    for(auto &r : mat) {
      if ((int)r.size() <= m) r.resize(m + 5);
    }
  }
};

// ============================================================================
//  Distance Implementations (Buffer-aware)
// ============================================================================

// 1. Levenshtein (Banded)
static inline int levenshtein_opt(const std::string &a, const std::string &b, int thr, std::vector<int>& row) {
  int n = a.size(), m = b.size();
  if (std::abs(n - m) > thr) return thr + 1;
  if ((int)row.size() <= m) row.resize(m + 1);
  
  std::iota(row.begin(), row.begin() + m + 1, 0);
  
  for (int i = 1; i <= n; ++i) {
    int start_j = std::max(1, i - thr);
    int end_j   = std::min(m, i + thr);

    // prev must hold the diagonal D[i-1][start_j-1]. When the band starts past
    // column 1 this is row[start_j-1] BEFORE it is overwritten by the boundary
    // guard - seeding it with D[i-1][0] (the old behaviour) corrupts the first
    // banded cell and makes the routine miss pairs whose only optimal alignment
    // rides the band edge (e.g. distance == |len diff| == thr).
    int prev = row[start_j - 1];
    row[0] = i;
    int row_min = i;

    // Fill valid band
    if (start_j > 1) row[start_j - 1] = thr + 2; // Guard boundary
    
    for (int j = start_j; j <= end_j; ++j) {
      int cur = row[j];
      int cost = (a[i-1] == b[j-1] ? 0 : 1);
      row[j] = std::min({ row[j-1] + 1, cur + 1, prev + cost });
      prev = cur;
      row_min = std::min(row_min, row[j]);
    }
    if (row_min > thr) return thr + 1;
  }
  return row[m];
}

// 3. Damerau-Levenshtein (Buffered)
static inline int damerau_opt(const std::string &a, const std::string &b, int thr, DPWorkspace &ws) {
  int n = a.size(), m = b.size();
  if (std::abs(n - m) > thr) return thr + 1;
  
  ws.resize_mat(n + 2, m + 2);
  auto &H = ws.mat;
  
  // Reset borders
  const int INF = thr + 2; // Optimization: Don't need full n+m
  H[0][0] = INF;
  for (int i = 0; i <= n + 1; ++i) { H[i][0] = INF; H[i][1] = i > 0 ? i-1 : 0; }
  for (int j = 0; j <= m + 1; ++j) { H[0][j] = INF; H[1][j] = j > 0 ? j-1 : 0; }
  
  // Char tracker (ASCII 0-255)
  if (ws.last_row_tracker.empty()) ws.last_row_tracker.resize(256, 0);
  std::fill(ws.last_row_tracker.begin(), ws.last_row_tracker.end(), 0);
  
  for (int i = 1; i <= n; ++i) {
    int last_match_col = 0;
    int row_min = INF;
    
    for (int j = 1; j <= m; ++j) {
      int i1 = i + 1; // 1-based index in H for string char i
      int j1 = j + 1;
      
      int last_match_row = ws.last_row_tracker[(unsigned char)b[j-1]];
      int cost = (a[i-1] == b[j-1] ? 0 : 1);
      
      int d_sub = H[i1-1][j1-1] + cost;
      int d_ins = H[i1][j1-1] + 1;
      int d_del = H[i1-1][j1] + 1;
      int d_trans = INF;
      
      if (last_match_row > 0 && last_match_col > 0) {
        d_trans = H[last_match_row][last_match_col] + (i - last_match_row - 1) + 1 + (j - last_match_col - 1);
      }
      
      H[i1][j1] = std::min({d_sub, d_ins, d_del, d_trans});
      if (cost == 0) last_match_col = j;
      row_min = std::min(row_min, H[i1][j1]);
    }
    ws.last_row_tracker[(unsigned char)a[i-1]] = i;
    if (row_min > thr) return thr + 1; 
  }
  return H[n+1][m+1];
}

// 4. Needleman-Wunsch (Buffered)
static inline int nw_opt(const std::string &a, const std::string &b, const SubstMatrix &mat, DPWorkspace &ws) {
  int n = a.size(), m = b.size();
  ws.resize_mat(n + 1, m + 1);
  auto &dp = ws.mat;
  
  // Init with proper gap penalties
  dp[0][0] = 0;
  for (int i = 1; i <= n; ++i) dp[i][0] = mat.gap_open + (i-1) * mat.gap_extend;
  for (int j = 1; j <= m; ++j) dp[0][j] = mat.gap_open + (j-1) * mat.gap_extend;
  
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int match = dp[i-1][j-1] + mat.get_score(a[i-1], b[j-1]);
      int del = dp[i-1][j] + mat.gap_extend; // Linear approximation for speed (vs affine state machine)
      int ins = dp[i][j-1] + mat.gap_extend;
      dp[i][j] = std::max({match, del, ins});
    }
  }
  // Convert to distance-like (inverted score)
  int max_s = std::max(n, m) * 5; 
  return std::max(0, max_s - dp[n][m]);
}

// 5. Smith-Waterman (Buffered)
static inline int sw_opt(const std::string &a, const std::string &b, const SubstMatrix &mat, DPWorkspace &ws) {
  int n = a.size(), m = b.size();
  ws.resize_mat(n + 1, m + 1);
  auto &dp = ws.mat;
  
  int max_val = 0;
  // Init (zeros for local alignment)
  for (int i = 0; i <= n; ++i) dp[i][0] = 0;
  for (int j = 0; j <= m; ++j) dp[0][j] = 0;
  
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int match = dp[i-1][j-1] + mat.get_score(a[i-1], b[j-1]);
      int del = dp[i-1][j] + mat.gap_extend;
      int ins = dp[i][j-1] + mat.gap_extend;
      dp[i][j] = std::max({0, match, del, ins});
      if (dp[i][j] > max_val) max_val = dp[i][j];
    }
  }
  int max_possible = std::min(n, m) * 5; 
  return std::max(0, max_possible - max_val);
}

// ============================================================================
//  Matrix Parser Helper
// ============================================================================
SubstMatrix parse_matrix(NumericMatrix mat, CharacterVector rnames, int go, int ge) {
  SubstMatrix sm;
  sm.gap_open = go;
  sm.gap_extend = ge;
  int n = mat.nrow();
  sm.scores.resize(n, std::vector<int>(n));
  
  for(int i=0; i<n; ++i) {
    std::string s = as<std::string>(rnames[i]);
    if(!s.empty()) sm.char_map[(unsigned char)s[0]] = i;
    for(int j=0; j<n; ++j) sm.scores[i][j] = (int)mat(i,j);
  }
  return sm;
}

// ============================================================================
//  Main Export
// ============================================================================
// [[Rcpp::export]]
DataFrame fast_edge_list(CharacterVector seqs,
                         double thresh = 1.0,
                         Nullable<CharacterVector> v_gene = R_NilValue,
                         Nullable<CharacterVector> j_gene = R_NilValue,
                         bool match_v = false,
                         bool match_j = false,
                         std::string metric = "levenshtein",
                         std::string normalize = "none",
                         Nullable<NumericMatrix> subst_matrix = R_NilValue,
                         int gap_open = -10,
                         int gap_extend = -1,
                         std::string expand = "clique")
{
  const bool star = (expand == "star");
  int n = seqs.size();
  if (n < 2) stop("Need at least 2 sequences.");
  
  // 1. Prepare Data
  std::vector<std::string> s(n);
  std::vector<int> lens(n);
  for(int i=0; i<n; ++i) { s[i] = as<std::string>(seqs[i]); lens[i] = s[i].size(); }
  
  std::vector<std::string> v_vec(n), j_vec(n);
  if(match_v && v_gene.isNotNull()) v_vec = as<std::vector<std::string>>(v_gene);
  if(match_j && j_gene.isNotNull()) j_vec = as<std::vector<std::string>>(j_gene);
  
  // 2. Deduplicate to representatives.
  //    A representative is a unique (sequence, V-used, J-used) tuple. V/J are
  //    only part of the key when the corresponding filter is on, so a rep never
  //    spans more than one (V,J) block. All distance work runs over reps; the
  //    member node indices are expanded back into edges afterwards.
  const char SEP = '\x1f';                 // unit separator: cannot occur in seq/gene
  std::unordered_map<std::string,int> rep_of;
  std::vector<std::vector<int>> rep_members;   // rep id -> node indices (ascending)
  std::vector<std::string>      rep_seq;
  std::vector<int>              rep_len;
  std::map<std::string, std::vector<int>> block_reps;  // block key -> rep ids

  rep_of.reserve(n * 2);
  for(int i=0; i<n; ++i) {
    std::string block_key;
    if (match_v) block_key += v_vec[i];
    block_key += SEP;
    if (match_j) block_key += j_vec[i];

    std::string full_key = s[i];
    full_key += SEP;
    full_key += block_key;

    auto it = rep_of.find(full_key);
    if (it == rep_of.end()) {
      int r = (int)rep_seq.size();
      rep_of.emplace(std::move(full_key), r);
      rep_seq.push_back(s[i]);
      rep_len.push_back(lens[i]);
      rep_members.push_back(std::vector<int>{i});
      block_reps[block_key].push_back(r);
    } else {
      rep_members[it->second].push_back(i);
    }
  }

  // 3. Prepare Substitution Matrix
  SubstMatrix smat;
  if (metric == "nw" || metric == "sw") {
    if (subst_matrix.isNotNull()) {
      NumericMatrix nm(subst_matrix);
      // FIXED: Ambiguous operator[] access
      List dimnames = nm.attr("dimnames");
      smat = parse_matrix(nm, dimnames[0], gap_open, gap_extend);
    } else {
      // Fallback Identity
      CharacterVector abc = CharacterVector::create("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
      NumericMatrix ident(20,20);
      
      // FIXED: Incorrect loop syntax/logic for matrix fill
      std::fill(ident.begin(), ident.end(), -1); // Fill all with mismatch (-1)
      for(int i=0; i<20; ++i) ident(i,i) = 1;    // Set diagonal to match (1)
      
      smat = parse_matrix(ident, abc, gap_open, gap_extend);
    }
  }
  
  // 4. Output Storage (0-based node indices; labels are attached in R)
  std::vector<int> out_i, out_k;
  std::vector<double> out_d;
  
  // 5. Per-pair decision (the exact, brute-force-faithful verifier).
  //    Encodes metric / normalize as ints so the hot path avoids string
  //    comparisons. Returns true and sets `fd_out` when an edge is kept,
  //    reproducing the original maxd / length-gate / normalize / post-gate
  //    logic byte-for-byte.
  enum Metric { LEV, HAM, DAM, MNW, MSW };
  enum Norm   { N_NONE, N_MAXLEN, N_MEANLEN };
  Metric mcode = (metric == "hamming") ? HAM :
                 (metric == "damerau") ? DAM :
                 (metric == "nw")      ? MNW :
                 (metric == "sw")      ? MSW : LEV;
  Norm   ncode = (normalize == "maxlen") ? N_MAXLEN :
                 (normalize == "length") ? N_MEANLEN : N_NONE;

  auto compute_pair = [&](const std::string& A, int la,
                          const std::string& B, int lb,
                          DPWorkspace& ws, double& fd_out) -> bool {
    int max_l = la > lb ? la : lb;
    if (max_l == 0) return false;                 // matches original `max_l==0 continue`
    int maxd;
    if (thresh >= 1.0) {
      maxd = (int)thresh;
    } else {
      double norm_len = (ncode == N_MEANLEN) ? (la + lb) / 2.0 : (double)max_l;
      maxd = (int)(thresh * norm_len + THRESH_EPS);
    }
    if (mcode != MNW && mcode != MSW && std::abs(la - lb) > maxd) return false;

    int d = 0;
    switch (mcode) {
      case LEV: d = levenshtein_opt(A, B, maxd, ws.v1); break;
      case HAM:
        if (la != lb) d = maxd + 1;
        else { for (int p = 0; p < la; ++p) if (A[p] != B[p]) { if (++d > maxd) break; } }
        break;
      case DAM: d = damerau_opt(A, B, maxd, ws); break;
      case MNW: d = nw_opt(A, B, smat, ws); break;
      case MSW: d = sw_opt(A, B, smat, ws); break;
    }
    if (d > maxd) return false;

    double fd = (double)d;
    if (ncode == N_MAXLEN)       fd /= max_l;
    else if (ncode == N_MEANLEN) fd /= ((la + lb) / 2.0);
    if (thresh < 1.0 && ncode != N_NONE && fd > thresh) return false;

    fd_out = fd;
    return true;
  };

  // 6. Per-block index: reps sorted by id, length buckets, and a block-level
  //    integer distance cap k_block. Because every per-pair maxd <= k_block
  //    (norm_len <= Lmax_block), restricting candidates to length-compatible
  //    buckets never drops a real edge - the exact per-pair verifier decides.
  int R = (int)rep_seq.size();
  struct BlockIdx {
    std::vector<int> rep_ids;                  // ascending rep id
    std::map<int, std::vector<int>> buckets;   // length -> ascending rep ids
    int kblock = 0;
    int Lmax = 0;
  };
  std::vector<BlockIdx> bidx;
  bidx.reserve(block_reps.size());
  std::vector<int> rep_block(R);
  {
    int b = 0;
    for (auto& kv : block_reps) {
      BlockIdx bi;
      bi.rep_ids = std::move(kv.second);
      int Lmax = 0;
      for (int r : bi.rep_ids) {
        rep_block[r] = b;
        bi.buckets[rep_len[r]].push_back(r);
        if (rep_len[r] > Lmax) Lmax = rep_len[r];
      }
      bi.kblock = (thresh >= 1.0) ? (int)thresh : (int)(thresh * (double)Lmax + THRESH_EPS);
      bi.Lmax = Lmax;
      bidx.push_back(std::move(bi));
      ++b;
    }
  }

  // 6b. Hamming pigeonhole index (one per equal-length bucket). Two equal-length
  //     sequences within Hamming distance k must share at least one of k+1
  //     positioned segments, so candidate partners are the reps that collide on
  //     a segment. Buckets where k >= L (every pair trivially passes) fall back
  //     to a plain bucket scan.
  struct SegBucket {
    bool brute = false;
    int nseg = 0;
    std::vector<int> seg_start, seg_len;
    std::unordered_map<std::string, std::vector<int>> index;  // pos-tag+seg -> ascending rep ids
  };
  std::vector<std::unordered_map<int, SegBucket>> ham(mcode == HAM ? bidx.size() : 0);
  if (mcode == HAM) {
    for (size_t b = 0; b < bidx.size(); ++b) {
      for (auto& lb : bidx[b].buckets) {
        int L = lb.first;
        const std::vector<int>& reps = lb.second;
        if (L == 0 || reps.size() < 2) continue;          // no possible edges
        int kseg = (thresh >= 1.0) ? (int)thresh : (int)(thresh * (double)L + THRESH_EPS);
        SegBucket sb;
        if (kseg >= L) { sb.brute = true; ham[b].emplace(L, std::move(sb)); continue; }
        sb.nseg = kseg + 1;
        int base = L / sb.nseg, rem = L % sb.nseg, pos = 0;
        for (int j = 0; j < sb.nseg; ++j) {
          int len = base + (j < rem ? 1 : 0);
          sb.seg_start.push_back(pos);
          sb.seg_len.push_back(len);
          pos += len;
        }
        sb.index.reserve(reps.size() * sb.nseg);
        for (int r : reps) {
          const std::string& sq = rep_seq[r];
          for (int j = 0; j < sb.nseg; ++j) {
            std::string key;
            key.reserve(sb.seg_len[j] + 1);
            key += (char)(j + 1);                          // position tag (letters are >= 'A', no clash)
            key.append(sq, sb.seg_start[j], sb.seg_len[j]);
            sb.index[key].push_back(r);                    // reps ascending -> postings ascending
          }
        }
        ham[b].emplace(L, std::move(sb));
      }
    }
  }

  // 6c. Deletion-neighborhood (SymSpell) index per block for edit-distance
  //     metrics with a small radius. Levenshtein uses radius = k_block; the
  //     guard keeps the variant explosion bounded, otherwise a block falls back
  //     to the length-blocked bucket scan. (Damerau, radius 2*k_block, is added
  //     alongside in step 8.)
  const int  K_CAP  = 3;
  const long MAXVAR = 4096;
  struct DelIdx {
    bool use = false;
    int radius = 0;
    std::unordered_map<std::string, std::vector<int>> index;  // variant -> ascending rep ids
  };
  std::vector<DelIdx> del(bidx.size());
  if (mcode == LEV || mcode == DAM) {
    std::unordered_set<std::string> vars;
    for (size_t b = 0; b < bidx.size(); ++b) {
      int k = bidx[b].kblock;
      if (k < 0 || k > K_CAP) continue;
      int radius = (mcode == DAM) ? 2 * k : k;   // transposition costs 2 in Levenshtein
      if (deletion_variant_bound(bidx[b].Lmax, radius) > MAXVAR) continue;
      DelIdx& di = del[b];
      di.use = true;
      di.radius = radius;
      di.index.reserve(bidx[b].rep_ids.size() * 4);
      for (int r : bidx[b].rep_ids) {
        gen_deletes(rep_seq[r], radius, vars);
        for (const std::string& v : vars) di.index[v].push_back(r);
      }
    }
  }

  // 7. Parallel over reps. Each rep emits its identical-member edges (clique or
  //    star) plus the edges to length-compatible reps with a higher id (so each
  //    unordered pair is generated exactly once). Expansion maps rep-pairs back
  //    to member node indices.
#pragma omp parallel
{
  DPWorkspace ws;
  std::vector<int> loc_f, loc_t;
  std::vector<double> loc_d;
  std::vector<int> cand;                  // reusable candidate-partner buffer
  std::unordered_set<std::string> vars;   // reusable deletion-variant buffer

  // verify a rep-pair and, if kept, expand to member node pairs
  auto verify_expand = [&](int r, int sv) {
    double fd;
    if (!compute_pair(rep_seq[r], rep_len[r], rep_seq[sv], rep_len[sv], ws, fd)) return;
    const std::vector<int>& mr = rep_members[r];
    const std::vector<int>& ms = rep_members[sv];
    if (star) {
      loc_f.push_back(mr[0]); loc_t.push_back(ms[0]); loc_d.push_back(fd);
    } else {
      for (int a : mr) for (int b : ms) {
        loc_f.push_back(a); loc_t.push_back(b); loc_d.push_back(fd);
      }
    }
  };
  // consider partners in an ascending bucket whose rep id exceeds r
  auto consider = [&](int r, const std::vector<int>& bucket) {
    for (auto it = std::upper_bound(bucket.begin(), bucket.end(), r);
         it != bucket.end(); ++it)
      verify_expand(r, *it);
  };

#pragma omp for schedule(dynamic, 32)
  for (int r = 0; r < R; ++r) {
    // intra-rep edges (identical sequences, mutual distance d_self).
    const std::vector<int>& mr = rep_members[r];
    double fd;
    if (mr.size() >= 2 &&
        compute_pair(rep_seq[r], rep_len[r], rep_seq[r], rep_len[r], ws, fd)) {
      if (star) {
        int hub = mr[0];
        for (size_t a = 1; a < mr.size(); ++a) {
          loc_f.push_back(hub); loc_t.push_back(mr[a]); loc_d.push_back(fd);
        }
      } else {
        for (size_t a = 0; a < mr.size(); ++a)
          for (size_t b = a + 1; b < mr.size(); ++b) {
            loc_f.push_back(mr[a]); loc_t.push_back(mr[b]); loc_d.push_back(fd);
          }
      }
    }

    // inter-rep candidates, restricted by the (V,J,length) index.
    const BlockIdx& bi = bidx[rep_block[r]];
    if (mcode == MNW || mcode == MSW) {
      consider(r, bi.rep_ids);                       // alignment scores: no length pruning
    } else if (mcode == HAM) {
      // Hamming: equal length only, candidates via pigeonhole segments.
      const std::unordered_map<int, SegBucket>& hmap = ham[rep_block[r]];
      auto hit = hmap.find(rep_len[r]);
      if (hit == hmap.end()) {
        // bucket had < 2 reps or empty strings -> no edges
      } else if (hit->second.brute) {
        auto bit = bi.buckets.find(rep_len[r]);
        if (bit != bi.buckets.end()) consider(r, bit->second);
      } else {
        const SegBucket& sb = hit->second;
        const std::string& sq = rep_seq[r];
        cand.clear();
        for (int j = 0; j < sb.nseg; ++j) {
          std::string key;
          key.reserve(sb.seg_len[j] + 1);
          key += (char)(j + 1);
          key.append(sq, sb.seg_start[j], sb.seg_len[j]);
          auto pit = sb.index.find(key);
          if (pit != sb.index.end())
            for (auto p = std::upper_bound(pit->second.begin(), pit->second.end(), r);
                 p != pit->second.end(); ++p)
              cand.push_back(*p);
        }
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
        for (int s : cand) verify_expand(r, s);
      }
    } else {
      // edit distance (Levenshtein / Damerau)
      const DelIdx& di = del[rep_block[r]];
      if (di.use) {
        // SymSpell: candidates are reps sharing a deletion variant
        gen_deletes(rep_seq[r], di.radius, vars);
        cand.clear();
        for (const std::string& v : vars) {
          auto it = di.index.find(v);
          if (it != di.index.end())
            for (auto p = std::upper_bound(it->second.begin(), it->second.end(), r);
                 p != it->second.end(); ++p)
              cand.push_back(*p);
        }
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
        for (int s : cand) verify_expand(r, s);
      } else {
        // fallback: length-blocked bucket scan, |dL| <= k_block
        int L = rep_len[r], k = bi.kblock;
        for (int Lp = L - k; Lp <= L + k; ++Lp) {
          auto it = bi.buckets.find(Lp);
          if (it != bi.buckets.end()) consider(r, it->second);
        }
      }
    }
  }
#pragma omp critical
{
  out_i.insert(out_i.end(), loc_f.begin(), loc_f.end());
  out_k.insert(out_k.end(), loc_t.begin(), loc_t.end());
  out_d.insert(out_d.end(), loc_d.begin(), loc_d.end());
}
}

  return DataFrame::create(_["i"] = out_i, _["k"] = out_k, _["dist"] = out_d, _["stringsAsFactors"] = false);
}
