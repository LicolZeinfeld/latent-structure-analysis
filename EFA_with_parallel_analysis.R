# ============================================================
# EFA_with_parallel_analysis.R
# 
# Purpose:
#   For humans and bots, separately:
#     1) Preprocess data
#     2) Compute tetrachoric correlation matrices
#     3) Retain factors using parallel analysis
#     4) Fit EFA using retained number of factors
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(psych)
})

# -----------------------------
# USER CONFIG
# -----------------------------
humans_csv <- ""
bots_csv   <- ""

rotate_method <- "oblimin"
fm_method     <- "minres"

# Parallel analysis settings
pa_n_iter <- 30

# Optional saving
save_plots <- TRUE
out_dir <- ""

# -----------------------------
# Helpers
# -----------------------------
read_binary_noheader <- function(path, prefix = "Q") {
  df <- read_csv(path, col_names = FALSE, show_col_types = FALSE)
  names(df) <- paste0(prefix, seq_len(ncol(df)))
  df <- as.data.frame(df)
  
  # force numeric 0/1 if possible
  df[] <- lapply(df, function(x) as.integer(as.character(x)))
  df
}

drop_zero_variance <- function(df) {
  keep <- sapply(df, function(x) var(x, na.rm = TRUE) > 0)
  df[, keep, drop = FALSE]
}

make_common_item_set <- function(df1, df2) {
  common <- intersect(names(df1), names(df2))
  list(
    df1 = df1[, common, drop = FALSE],
    df2 = df2[, common, drop = FALSE]
  )
}

get_tetra_R <- function(df) {
  psych::tetrachoric(df)$rho
}

# -----------------------------
# load data
# -----------------------------
hum <- read_binary_noheader(humans_csv, prefix = "Q")
bot <- read_binary_noheader(bots_csv, prefix = "Q")

# manually drop selected items from instrument and exclude from analyses
items_to_drop <- c()

hum <- hum[, !(names(hum) %in% items_to_drop), drop = FALSE]
bot <- bot[, !(names(bot) %in% items_to_drop), drop = FALSE]

# remove zero-variance items separately
hum <- drop_zero_variance(hum)
bot <- drop_zero_variance(bot)

# keep only common items across both groups
aligned <- make_common_item_set(hum, bot)
hum <- aligned$df1
bot <- aligned$df2

cat("Final #items used (common, nonzero variance):", ncol(hum), "\n\n")

if (ncol(hum) < 3) stop("Too few usable items after filtering.")

# -----------------------------
# tetrachoric correlations
# -----------------------------
R_h <- get_tetra_R(hum)
R_b <- get_tetra_R(bot)

# -----------------------------
# parallel analysis factor retention
# -----------------------------
pa_h <- psych::fa.parallel(
  hum,
  fm = fm_method,
  fa = "fa",
  cor = "tet",
  n.iter = pa_n_iter,
  plot = FALSE
)

pa_b <- psych::fa.parallel(
  bot,
  fm = fm_method,
  fa = "fa",
  cor = "tet",
  n.iter = pa_n_iter,
  plot = FALSE
)

K_h <- pa_h$nfact
K_b <- pa_b$nfact

if (is.null(K_h) || is.na(K_h) || K_h < 1) K_h <- 1
if (is.null(K_b) || is.na(K_b) || K_b < 1) K_b <- 1

cat("Parallel analysis suggested K (humans):", K_h, "\n")
cat("Parallel analysis suggested K (bots):  ", K_b, "\n\n")

# -----------------------------
# fit EFA for humans and bots exactly at the selected K from parallel analysis
# -----------------------------
efa_h <- psych::fa(
  r = R_h,
  nfactors = K_h,
  n.obs = nrow(hum),
  rotate = rotate_method,
  fm = fm_method
)

efa_b <- psych::fa(
  r = R_b,
  nfactors = K_b,
  n.obs = nrow(bot),
  rotate = rotate_method,
  fm = fm_method
)

# -----------------------------
# optional output directory
# -----------------------------
if (save_plots) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

