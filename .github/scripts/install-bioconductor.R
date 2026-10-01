# Run after the standard package updater, which may leave Bioconductor out.
Sys.setenv(RENV_PROFILE = "lesson-requirements")
source("renv/activate.R")
lesson_library <- renv::paths$library()
.libPaths(c(lesson_library, .libPaths()))
options(timeout = 1200, warn = 1)

required <- c(
  "AnnotationHub", "clusterProfiler", "ComplexHeatmap", "DMRcate", "GEOquery",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "IlluminaHumanMethylationEPICv2manifest", "limma", "methylKit", "minfi",
  "missMethyl", "org.Hs.eg.db"
)

renv::install("BiocManager", library = lesson_library, prompt = FALSE)
version <- renv::settings$bioconductor.version()
if (is.null(version) || !nzchar(version)) {
  version <- as.character(BiocManager::version())
}
renv::settings$bioconductor.version(version)
repos <- BiocManager::repositories(version = version)
options(repos = repos)
cat("Installing Bioconductor", version, "into", lesson_library, "\n")
print(repos)

# Explicit sources avoid relying on the updater's automatic package discovery.
renv::install(paste0("bioc::", required), library = lesson_library,
  repos = repos, prompt = FALSE)

loadable <- vapply(required, requireNamespace, logical(1), quietly = TRUE,
  lib.loc = lesson_library)
if (!all(loadable)) {
  stop("Bioconductor packages failed to load: ",
    paste(required[!loadable], collapse = ", "))
}

lockfile <- renv::paths$lockfile()
previous <- renv::lockfile_read(lockfile)
# Preserve the updater's package records and include all recursive dependencies.
# Include R's recommended packages from the active system/sandbox library too.
# Restricting this to lesson_library incorrectly reports them as uninstalled.
renv::snapshot(library = .libPaths(), lockfile = lockfile,
  packages = unique(c(names(previous$Packages), "BiocManager", required)),
  repos = repos, prompt = FALSE)
lock <- renv::lockfile_read(lockfile)
missing <- setdiff(required, names(lock$Packages))
if (length(missing)) {
  stop("Lockfile is missing: ", paste(missing, collapse = ", "))
}
if (is.null(lock$Bioconductor$Version)) {
  stop("Lockfile does not record a Bioconductor version")
}
cat("All", length(required), "Bioconductor packages load and are recorded.\n")
