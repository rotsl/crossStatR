# Render the public reports alongside the pkgdown site.
# Usage: Rscript scripts/build_reports.R <site-directory>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[1L])) {
  stop("An existing site directory is required.")
}
destination <- file.path(normalizePath(args[1L]), "docs")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
reports <- c("BENCHMARK_REPORT", "VALIDATION_REPORT", "METHODS", "REFERENCES")
for (name in reports) {
  output <- rmarkdown::render(
    file.path("docs", paste0(name, ".md")),
    output_format = rmarkdown::html_document(self_contained = TRUE, toc = TRUE),
    output_file = paste0(name, ".html"), output_dir = destination,
    intermediates_dir = tempdir(), quiet = TRUE)
  html <- readLines(output, warn = FALSE)
  html <- gsub('href="([^":]+)[.]md([#][^"]*)?"',
               'href="\\1.html\\2"', html, perl = TRUE)
  writeLines(html, output)
}
results <- list.files("docs", pattern = "[.]csv$", full.names = TRUE)
stopifnot(all(file.copy(results, destination, overwrite = TRUE)),
          all(file.exists(file.path(destination, paste0(reports, ".html")))))
