# ---------------------------------------------
# Simple MOI -> infection count calculator
# Assumption: phage adsorptions per cell ~ Poisson(lambda = MOI)
# ---------------------------------------------

# ---- Inputs you change ----
N_cells <- 1e9        # total number of cells
moi_values <- seq(0, 5, by = 0.5)  # MOIs to evaluate
k_exact <- 1          # "infected by exactly k phages"
x_or_more <- 3        # "infected by x or more phages"
K_max <- 15           # maximum k shown in histogram

# ---- Helper functions ----

# Expected cells with exactly k adsorptions
cells_exact_k <- function(moi, N, k) {
  N * dpois(k, lambda = moi)
}

# Expected cells not infected (k = 0)
cells_not_infected <- function(moi, N) {
  N * dpois(0, lambda = moi)
}

# Expected cells with >= x adsorptions
cells_x_or_more <- function(moi, N, x) {
  N * ppois(x - 1, lambda = moi, lower.tail = FALSE)
}

# ---- Calculate across MOIs ----
results <- data.frame(
  MOI = moi_values,
  not_infected = cells_not_infected(moi_values, N_cells),
  infected_exact_k = cells_exact_k(moi_values, N_cells, k_exact),
  infected_x_or_more = cells_x_or_more(moi_values, N_cells, x_or_more)
)

print(results)

# ---------------------------------------------
# Histogram for one selected MOI
# ---------------------------------------------

# Choose which MOI to visualize
moi_plot <- 1

# Compute expected counts for k = 0..K_max
k_vals <- 0:K_max
expected_counts <- N_cells * dpois(k_vals, lambda = moi_plot)

# Add final bin for >= K_max+1
expected_counts <- c(expected_counts,
                     N_cells * ppois(K_max, lambda = moi_plot, lower.tail = FALSE))
k_labels <- c(k_vals, paste0(">=", K_max + 1))

# ---- Plot as discrete histogram (barplot) ----
barplot(expected_counts,
        names.arg = k_labels,
        col = "steelblue",
        border = "white",
        xlab = "Number of phages adsorbed per cell (k)",
        ylab = "Expected number of cells",
        main = paste("Distribution of Phage Adsorption\nMOI =", moi_plot))

######
# ---------------------------------------------
# Overlapping MOI distributions
# Assumption: adsorptions per cell ~ Poisson(lambda = MOI)
# ---------------------------------------------

# ---- Inputs you change ----
N_cells <- 1e8                    # total number of cells
moi_values <- seq(0.5, 5, by=1) # MOIs to compare
K_max <- 10                       # maximum k shown on x-axis

# ---- Prepare x-axis (number of adsorbed phages) ----
k_vals <- 0:K_max

# ---- Set up empty plot ----
plot(NULL,
     xlim = c(0, K_max),
     ylim = c(0, 8e07),
     xlab = "Number of phages adsorbed per cell (k)",
     ylab = "Expected number of cells",
     main = "Overlapping Phage Adsorption Distributions")

# ---- Add one line per MOI ----
colors <- rainbow(length(moi_values))  # distinct colors

for (i in seq_along(moi_values)) {
  
  moi <- moi_values[i]
  
  # Expected counts for this MOI
  expected_counts <- N_cells * dpois(k_vals, lambda = moi)
  
  # Add line
  lines(k_vals,
        expected_counts,
        type = "l",
        lwd = 2,
        col = colors[i])
}

# ---- Add legend ----
legend("topright",
       legend = paste("MOI =", moi_values),
       col = colors,
       lwd = 2,
       cex = 0.8)
########

# -------------------------------------------------------
# "Majority have >1 phage adsorbed" threshold diagram
# Assumption: adsorptions per cell K ~ Poisson(lambda = MOI)
# We want the MOI where P(K >= 2) = 0.5
# -------------------------------------------------------

# ---- Fraction function: P(K >= 2) ----
p_ge_2 <- function(moi) {
  1 - ppois(1, lambda = moi)  # 1 - P(K <= 1)
}

# ---- Solve for MOI where P(K >= 2) = 0.5 ----
# We solve p_ge_2(moi) - 0.5 = 0 over a reasonable interval
threshold <- uniroot(function(m) p_ge_2(m) - 0.5, interval = c(0, 10))$root

# ---- Make a curve over MOI values ----
moi_grid <- seq(0, 6, by = 0.01)

p0   <- dpois(0, lambda = moi_grid)           # P(K = 0)
p1   <- dpois(1, lambda = moi_grid)           # P(K = 1)
pge2 <- 1 - ppois(1, lambda = moi_grid)       # P(K >= 2)

# ---- Plot the diagram (base R) ----
plot(moi_grid, pge2, type = "l", lwd = 2,
     xlab = "MOI (mean phage adsorptions per cell)",
     ylab = "Fraction of cells",
     ylim = c(0, 1),
     main = "When does a majority of cells have >1 phage adsorbed?")

# Add helpful reference curves (optional but informative)
lines(moi_grid, p0, lwd = 2, lty = 2)  # P(K=0)
lines(moi_grid, p1, lwd = 2, lty = 3)  # P(K=1)

# Horizontal line at 0.5 (majority threshold)
abline(h = 0.5, lty = 2)

# Vertical line at the threshold MOI
abline(v = threshold, lwd = 2, lty = 1)

# Mark the threshold point on the P(K>=2) curve
points(threshold, 0.5, pch = 19)

# Label it
text(threshold, 0.55,
     labels = paste0("P(K>=2)=0.5 at MOI ≈ ", round(threshold, 3)),
     pos = 4)

legend("bottomright",
       legend = c("P(K >= 2)", "P(K = 0)", "P(K = 1)"),
       lty = c(1, 2, 3),
       lwd = 2,
       bty = "n")

# Print the threshold value to the console too
cat("MOI threshold where majority have >1 adsorption (P(K>=2)=0.5):",
    threshold, "\n")

