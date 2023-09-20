lib    <- "/tmp/rpkgs/build"
mirror <- "http://cran.r-project.org"
Ncpus  <- 8

# Install remotes, so we can pin dependency versions
install.packages("remotes", repos=mirror, Ncpus=Ncpus)
library(remotes)

# Install dependencies
# NOTE While this works, it seems to do a lot of unnecessary work as
# common dependencies are rebuilt for some reason. It might be worth
# listing all dependencies and transitive dependencies and topologically
# sorting them, to avoid this.
requirements <- rbind(
  c("codetools",  "0.2-20"),
  c("dplyr",      "1.1.4"),
  c("pheatmap",   "1.0.12"),
  c("rmarkdown",  "2.29"),
  c("tidyr",      "1.3.1"),
  c("timeSeries", "4041.111"),
  c("ggplot2",    "3.5.1")
)

requirements <- as.data.frame(requirements)
colnames(requirements) <- c("package", "version")

apply(
  requirements,
  1,
  function(row) {
    install_version(
      package=row["package"],
      version=row["version"],
      lib=lib,
      repos=mirror,
      Ncpus=Ncpus
    )
  }
)
