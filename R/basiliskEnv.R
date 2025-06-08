#' @import basilisk
immApexEnv <- BasiliskEnvironment(
  envname = "immApexEnv",
  pkgname = "immApex",
  packages = c(
    "python=3.10",
    "pip=23.3.1",
    "setuptools=69.2.*",
    "wheel=0.42.*",
    "numpy=1.26", 
    "h5py=3.10", 
    "keras=3.3.*",
    "tensorflow=2.16.*"
  )
)