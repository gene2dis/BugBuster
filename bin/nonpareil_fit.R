#!/usr/bin/env Rscript
# Nonpareil curve fit of one sample (NONPAREIL module): Nonpareil.curve() from the Nonpareil R package on the
# .npo the nonpareil binary wrote, and two tables for the run summary, so that no other script re-implements the
# coverage correction or the model.
#
# Usage: nonpareil_fit.R <sample>.npo <sample> [reads_shorter_than_k]
#
# <sample>_nonpareil.tsv (one row):
#   sample, reads (R: reads Nonpareil used; R1 only, so read pairs), reads_shorter_than_k (R1 reads the module
#   left out, shorter than the k-mer; blank when not given), read_length (L), al (AL: the read length the
#   package scales the effort by), kappa (redundancy at the sequenced depth), coverage (C: estimated average coverage
#   at the sequenced depth), diversity (Nd), lr_star (LR*: effort for 95 % coverage on the fit's axis, AL x reads),
#   depth_95pct_pairs (LR* / AL: read pairs for 95 % coverage), model_r (correlation of the model with the curve),
#   has_model,
#   consistent (the package's curve checks), warning (the package's warnings, '; '-separated)
# <sample>_nonpareil_curve.tsv: the curve as Nonpareil.curve reads it (kind 'observed': x_adj, coverage, its
#   25th-75th percentiles) and the model's prediction on a log grid (kind 'model'), x also in read pairs
#   (x_adj / AL, so the observed curve ends at R pairs).
#
# Nd, LR* and the model points are blank when the curve fails the package's checks (e.g. too few reads) or the
# model does not converge: the coverage is still reported. Nd is also blank when the fitted shape a <= 1 (the
# package then leaves it at 0, which is not a diversity value).
#
# AL: the package reads a k-mer .npo's ksize as text and tests it against 1001 as text ("24" < "1001" is false),
# so with k = 24 it treats the file as an alignment-kernel curve and AL = L instead of L - k + 1 (package 3.5.3).
# The fit's axis is then L x reads, the same as LR = R x L. Dividing by the package's own AL gives read pairs
# whichever way it reads the file.

suppressPackageStartupMessages(library(Nonpareil))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 2:3) stop("usage: nonpareil_fit.R <sample>.npo <sample> [reads_shorter_than_k]")
npo <- args[1]
sample <- args[2]
short_reads <- if (length(args) == 3) args[3] else ""

warns <- character(0)
np <- withCallingHandlers(
    Nonpareil.curve(npo, plot = FALSE, label = sample, col = "black"),
    warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
)

num <- function(x) if (length(x) == 1 && is.finite(x)) x else NA_real_
fmt <- function(x) ifelse(is.na(x), "", format(x, digits = 10, scientific = FALSE, trim = TRUE))

has_model <- isTRUE(np$has.model)
if (has_model && coef(np$model)[["a"]] <= 1) {
    warns <- c(warns, "Fitted shape a <= 1: no diversity index.")
}
diversity <- if (has_model && coef(np$model)[["a"]] > 1) num(np$diversity) else NA_real_
lr_star <- if (has_model) num(np$LRstar) else NA_real_
al <- num(np$AL)

values <- data.frame(
    sample = sample,
    reads = fmt(num(np$R)),
    reads_shorter_than_k = short_reads,
    read_length = fmt(num(np$L)),
    al = fmt(al),
    kappa = fmt(num(np$kappa)),
    coverage = fmt(num(np$C)),
    diversity = fmt(diversity),
    lr_star = fmt(lr_star),
    depth_95pct_pairs = fmt(if (is.na(lr_star) || is.na(al)) NA_real_ else lr_star / al),
    model_r = fmt(if (has_model) num(np$modelR) else NA_real_),
    has_model = if (has_model) "true" else "false",
    consistent = if (isTRUE(np$consistent)) "true" else "false",
    warning = paste(unique(gsub("[\t\r\n]+", " ", warns)), collapse = "; "),
    check.names = FALSE
)
write.table(values, paste0(sample, "_nonpareil.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

curve <- data.frame(
    sample = sample, kind = "observed", x_adj = np$x.adj, x_pairs = np$x.adj / al,
    coverage = np$y.cov, coverage_p25 = np$y.p25, coverage_p75 = np$y.p75
)
if (has_model) {
    # out to past the sequenced depth and LR*, so the figure shows the projection to 95 % coverage
    lo <- min(np$x.adj[np$x.adj > 0])
    hi <- max(10 * max(np$x.adj), if (is.na(lr_star)) 0 else 1.5 * lr_star)
    grid <- exp(seq(log(lo), log(hi), length.out = 200))
    curve <- rbind(curve, data.frame(
        sample = sample, kind = "model", x_adj = grid, x_pairs = grid / al,
        coverage = as.numeric(predict(np, lr = grid)), coverage_p25 = NA_real_, coverage_p75 = NA_real_
    ))
}
for (col in c("x_adj", "x_pairs", "coverage", "coverage_p25", "coverage_p75")) curve[[col]] <- fmt(curve[[col]])
write.table(curve, paste0(sample, "_nonpareil_curve.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (length(warns)) message("Nonpareil.curve: ", paste(unique(warns), collapse = "; "))
