# ============================================================
# phiKZ_fitness_histogram.R
#
# Reproduces the fitness distribution histogram from Fig S3
# of Chan et al. 2025 (bioRxiv 10.1101/2025.11.23.690004)
#
# Input:  s3_curve_corrected.json  (pixel-digitised rank curve)
#         OR your own data frame with columns:
#           gene_id, insertion_index  (linear scale)
#           -- or --
#           gene_id, log10_insertion_index
#
# Output: histogram showing the bimodal fitness distribution
# ============================================================

library(jsonlite)
library(ggplot2)
library(dplyr)

# ----------------------------------------------------------
# 1. LOAD DATA
#    If you have the raw TraDIS/TnSeq output instead, replace
#    this block with read.csv() / read.table() as appropriate
# ----------------------------------------------------------

raw <- fromJSON("s3_curve_corrected.json")   # data frame with gene_rank & log10_insertion_index

# ----------------------------------------------------------
# 2. CLEAN THE DIGITISED CURVE
#
#    Two artefacts from pixel extraction to remove:
#      (a) ranks 0-16: top-border line of the figure → log10 ≈ -0.94
#      (b) the essential-gene floor is plotted at log10 < -3.0
#          (below the visible axis); we keep those as-is
# ----------------------------------------------------------

clean <- raw %>%
  filter(log10_insertion_index < -0.95,   # remove top-border artefact
         gene_rank >= 16)                  # same artefact cuts off here

# Re-zero and rescale ranks to [0, 371]
clean <- clean %>%
  mutate(rank_scaled = (gene_rank - min(gene_rank)) /
                       (max(gene_rank) - min(gene_rank)) * 371)

# ----------------------------------------------------------
# 3. INTERPOLATE TO ONE VALUE PER GENE (n = 371)
#
#    The digitised curve has ~600 x-samples for 371 genes.
#    We split into:
#      • essential pool  – points below the axis floor (< -3.0)
#      • dispensable set – the rising portion (≥ -3.0)
#    and interpolate the rising portion to exactly n_dispensable
#    evenly spaced gene ranks.
# ----------------------------------------------------------

# Fraction of curve at floor → estimate of essential gene count
n_total      <- 371
frac_ess     <- mean(clean$log10_insertion_index < -3.0)
n_essential  <- round(frac_ess * n_total)          # ≈ 107
n_dispensable <- n_total - n_essential

cat(sprintf("Essential genes (at floor): %d\n", n_essential))
cat(sprintf("Dispensable genes:          %d\n", n_dispensable))

# Interpolate rising portion
rising <- clean %>%
  filter(log10_insertion_index >= -3.0) %>%
  arrange(rank_scaled)

disp_ranks <- seq(min(rising$rank_scaled), max(rising$rank_scaled),
                  length.out = n_dispensable)

disp_log10 <- approx(rising$rank_scaled,
                     rising$log10_insertion_index,
                     xout = disp_ranks)$y

disp_log10 <- pmin(pmax(disp_log10, -3.0), -0.75)   # clamp to axis range

# ----------------------------------------------------------
# 4. BUILD PLOTTING DATA FRAMES
# ----------------------------------------------------------

# Dispensable genes (histogram)
disp_df <- data.frame(log10_ii = disp_log10)

# Essential genes represented as a single bar to the left of the axis
ess_bar <- data.frame(
  xmin = -3.35, xmax = -3.13,
  ymin = 0,     ymax = n_essential,
  label = sprintf("Essential\n(n ≈ %d)", n_essential)
)

# Reference lines
vlines <- data.frame(
  xint  = c(-3.0, -2.0),
  label = c("Axis floor  (log₁₀ = −3)", "log₁₀ = −2  (ins. idx = 0.01)"),
  col   = c("#8B0000", "darkorange")
)

# ----------------------------------------------------------
# 5. PLOT
# ----------------------------------------------------------

p <- ggplot(disp_df, aes(x = log10_ii)) +

  # --- main histogram (dispensable genes) ---
  geom_histogram(
    aes(fill = after_stat(x)),
    binwidth = 0.12,
    boundary = -3.0,
    color    = "white",
    linewidth = 0.3
  ) +
  scale_fill_viridis_c(
    option   = "plasma",
    direction = -1,
    limits   = c(-3.3, -0.75),
    guide    = "none"
  ) +

  # --- essential-gene bar (off-axis, shown in dark red) ---
  geom_rect(
    data = ess_bar,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill  = "#8B0000",
    color = "white",
    linewidth = 0.4
  ) +
  annotate("text",
    x     = mean(c(ess_bar$xmin, ess_bar$xmax)),
    y     = n_essential + 2,
    label = ess_bar$label,
    size  = 3, color = "#8B0000", fontface = "bold", hjust = 0.5
  ) +

  # --- reference vertical lines ---
  geom_vline(xintercept = -3.0,  linetype = "dashed",
             color = "#8B0000",   linewidth = 0.8, alpha = 0.8) +
  geom_vline(xintercept = -2.0,  linetype = "dashed",
             color = "darkorange", linewidth = 0.8, alpha = 0.8) +

  # --- zone annotations ---
  annotate("text", x = -1.55, y = 40,
           label = "Dispensable\n(many insertions)",
           size = 3, color = "steelblue", hjust = 0.5) +
  annotate("text", x = -2.5, y = 8,
           label = "Moderate\nfitness effect",
           size = 3, color = "darkorange4", hjust = 0.5) +

  # --- axes and theme ---
  scale_x_continuous(
    name   = expression(log[10](insertion~index)),
    limits = c(-3.45, -0.60),
    breaks = c(-3.0, -2.5, -2.0, -1.5, -1.0),
    labels = function(x) sprintf("%.1f\n(%.3f)", x, 10^x)
  ) +
  scale_y_continuous(name = "Number of genes", expand = expansion(mult = c(0, 0.08))) +

  labs(
    title    = expression(paste(phi, "KZ gene fitness distribution")),
    subtitle = "Derived from Fig S3 rank-order curve  •  n = 371 genes",
    caption  = "Chan et al. 2025, bioRxiv 10.1101/2025.11.23.690004"
  ) +

  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 10),
    plot.caption  = element_text(color = "grey55", size = 8),
    axis.line     = element_line(color = "grey30"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

print(p)

# Save
ggsave("phiKZ_fitness_histogram.pdf", p, width = 8, height = 5.5)
ggsave("phiKZ_fitness_histogram.png", p, width = 8, height = 5.5, dpi = 200)
message("Saved phiKZ_fitness_histogram.pdf / .png")

# ----------------------------------------------------------
# 6. SUMMARY TABLE
# ----------------------------------------------------------

breaks <- c(-Inf, -3.0, -2.5, -2.0, -1.5, -1.0, Inf)
labels <- c("Essential (floor)", "–3.0 to –2.5", "–2.5 to –2.0",
            "–2.0 to –1.5", "–1.5 to –1.0", "> –1.0")

all_vals <- c(rep(-3.5, n_essential), disp_log10)   # -3.5 as floor placeholder

summary_tbl <- data.frame(log10_ii = all_vals) %>%
  mutate(category = cut(log10_ii, breaks = breaks, labels = labels,
                        right = FALSE, include.lowest = TRUE)) %>%
  count(category, name = "n_genes") %>%
  mutate(pct = round(100 * n_genes / sum(n_genes), 1))

print(summary_tbl)
