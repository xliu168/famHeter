#!/usr/bin/env Rscript
# Save the software versions used for a reproducible run.
outfile <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  "sessionInfo.txt"
}
capture.output(sessionInfo(), file = outfile)
cat("Session information saved to:", outfile, "\n")
