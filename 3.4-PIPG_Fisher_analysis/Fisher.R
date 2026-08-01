library(dplyr)
library(readr)

# -------------------------
# Fisher test helper
# -------------------------
fisher_pairwise <- function(df1, df2, gene_col, label, group1, group2) {
  
  a <- sum(df1[[gene_col]])
  b <- sum(df1$ORF_Count) - a
  c <- sum(df2[[gene_col]])
  d <- sum(df2$ORF_Count) - c
  
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  rownames(mat) <- c(group1, group2)
  colnames(mat) <- c("Functional", "Non_functional")
  
  ft <- fisher.test(mat)
  
  data.frame(
    Category = label,
    Comparison = paste(group1, "vs", group2),
    Functional_1 = a,
    Functional_2 = c,
    Odds_Ratio = ft$estimate,
    P_value = ft$p.value
  )
}

# -------------------------
# Load data
# -------------------------
pipg_vf  <- read_tsv("plasmid_30_new_vf.txt", show_col_types = FALSE)
rand_vf  <- read_tsv("plasmid_random_1000_vf.txt", show_col_types = FALSE)

pipg_arg <- read_tsv("plasmid_30_new_arg.txt", show_col_types = FALSE)
rand_arg <- read_tsv("plasmid_random_1000_arg.txt", show_col_types = FALSE)

pipg_mge <- read_tsv("plasmid_30_new_te.txt", show_col_types = FALSE)
rand_mge <- read_tsv("plasmid_random_1000_mge.txt", show_col_types = FALSE)


# -------------------------
# Split random plasmids by size
# -------------------------
split_random <- function(df) {
  list(
    all   = df,
    mega  = df %>% filter(Length > 100000),
    small = df %>% filter(Length <= 100000)
  )
}

rand_vf_s  <- split_random(rand_vf)
rand_arg_s <- split_random(rand_arg)
rand_te_s <- split_random(rand_te)

# -------------------------
# Run Fisher tests
# -------------------------
results <- bind_rows(
  
  # VF
  fisher_pairwise(pipg_vf, rand_vf_s$all,   "VF_Count",  "VF",  "PIPGs", "Random_all"),
  fisher_pairwise(pipg_vf, rand_vf_s$mega,  "VF_Count",  "VF",  "PIPGs", "Random_mega"),
  fisher_pairwise(pipg_vf, rand_vf_s$small, "VF_Count",  "VF",  "PIPGs", "Random_small"),
  
  # ARG
  fisher_pairwise(pipg_arg, rand_arg_s$all,   "ARG_Count", "ARG", "PIPGs", "Random_all"),
  fisher_pairwise(pipg_arg, rand_arg_s$mega,  "ARG_Count", "ARG", "PIPGs", "Random_mega"),
  fisher_pairwise(pipg_arg, rand_arg_s$small, "ARG_Count", "ARG", "PIPGs", "Random_small"),
  
  # TE
  fisher_pairwise(pipg_mge, rand_mge_s$all,   "MGE_Count", "TE", "PIPGs", "Random_all"),
  fisher_pairwise(pipg_mge, rand_mge_s$mega,  "MGE_Count", "TE", "PIPGs", "Random_mega"),
  fisher_pairwise(pipg_mge, rand_mge_s$small, "MGE_Count", "TE", "PIPGs", "Random_small"),

)

# -------------------------
# Multiple testing correction
# -------------------------
results$P_adj_BH <- p.adjust(results$P_value, method = "BH")

# -------------------------
# Save results
# -------------------------
write_tsv(results, "Fisher_exact_PIPG_random_all_and_size.tsv")

print(results)
