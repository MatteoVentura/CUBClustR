# ==============================================================================
# MLCCUB Analysis: University Satisfaction Dataset
# ==============================================================================
# This script demonstrates the use of the MLCCUB package to analyze the 
# 'univer' dataset from the CUB package. We investigate latent structures 
# in student satisfaction across 5 items (7-point scale).

# ------------------------------------------------------------------------------
# 0. SETUP AND LIBRARIES
# ------------------------------------------------------------------------------
# Please extract the .zip file to get the .tar.gz installation file
# install.packages("path/to/CUBClustR.tar.gz", repos = NULL, type = "source")
library(CUBClustR)
library(CUB)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. DATA PREPARATION
# ------------------------------------------------------------------------------
# Loading the 'univer' dataset: Evaluation of the University of Naples Federico II
data(univer)

# Selection of 5 relevant ordinal items (7-point scale)
# Items: Global Satisfaction, Willingness to enroll again, etc.
R <- univer[, 8:12]
R <- apply(R, 2, as.integer)

# Maximum categories for each item (m = 7)
m_vector <- rep(7, ncol(R))

# ------------------------------------------------------------------------------
# 2. MODEL SELECTION (BIC)
# ------------------------------------------------------------------------------
# We test from 1 to 7 latent classes to find the optimal balance 
# between model complexity and fit.
model_selection <- compare_MLCCUB_models(
  data     = R, 
  m        = m_vector, 
  K_vector = 1:7, 
  EM_iter  = 5, 
  max_iter = 1000
)

# Plotting the BIC trend
ggplot(model_selection$summary, aes(x = K, y = BIC)) +
  geom_line(color = "#2C3E50", linewidth = 1.2) +
  geom_point(color = "#BD4F6C", size = 3) +
  theme_bw() +
  labs(title = "Model Selection: BIC Trend (University Data)",
       x     = "Number of Latent Classes (K)",
       y     = "BIC") +
  scale_x_continuous(breaks = 1:7)

# ------------------------------------------------------------------------------
# 3. STABILITY & IDENTIFIABILITY CHECK
# ------------------------------------------------------------------------------
# Based on BIC results, we evaluate the 
# identifiability of a 5-cluster solution (K=5).
identifiability <- check_identifiability(
  data      = R, 
  m         = m_vector, 
  K         = 5, 
  n_boot    = 100, 
  n_cores   = 2,
  max_iter  = 1000, 
  EM_iter   = 5
)

# Display the ARI distribution: values near 1 indicate highly stable clusters
identifiability$plot

# ------------------------------------------------------------------------------
# 4. FINAL MODEL VISUALIZATION (K = 3)
# ------------------------------------------------------------------------------
# Extracting the 5-cluster model for interpretation
model_k5 <- model_selection$all_models$K3

# Visualization 1: Parameter Space (Feeling vs. Uncertainty)
# Allows identifying "satisfied", "unhappy", or "uncertain" groups
plot_params <- plot_MLCCUB_clusters(model_k5)
plot_params + labs(title = "University Clusters: Parameter Space")

# Visualization 2: Observed vs Theoretical Fit
# Checking how well the model mimics the real data distributions
plot_fit <- plot_CUB_fit(
  model_obj  = model_k5, 
  data       = R, 
  m          = m_vector,
  line_color = "#93B5C6"
)
plot_fit

# ------------------------------------------------------------------------------
# 5. GOODNESS OF FIT DIAGNOSTICS
# ------------------------------------------------------------------------------
# Quantitative assessment of the fit quality
fit_metrics <- calculate_WAD(model_obj = model_k5, data = R, m = m_vector)

# Diagnostic results
fit_metrics$WAD                   # Global index (closer to 0 is better)
fit_metrics$Cluster_Dissimilarity # Which cluster is harder to fit?
fit_metrics$Item_Dissimilarity    # Which item is less explained by the model?


# ------------------------------------------------------------------------------
# 6. NON-METRIC MULTIDIMENSIONAL SCALING (NMDS) - MANHATTAN
# ------------------------------------------------------------------------------
# We use NMDS to project respondents into a 2D space. 
# Unlike classical MDS, NMDS preserves the rank-order of distances, 
# making it ideal for non-linear ordinal rating data.

library(vegan) # Required for metaMDS

# 1. Run NMDS using Manhattan distance on a subsample
# k = 2 (dimensions), trymax = 50 (iterations to find best solution)
set.seed(123) # For reproducibility of the NMDS starts
idx_sample <- sample(1:nrow(R), 800) 
R_sub <- R[idx_sample, ]
cluster_sub <- model_k5$class[idx_sample]

nmds_res <- metaMDS(R_sub, 
                    distance = "manhattan", 
                    k = 2, 
                    trymax = 20,
                    autotransform = FALSE, 
                    trace = 1)

# 2. Check the Stress value (A measure of how well the 2D map fits the data)
# Stress < 0.1: Great fit; < 0.2: Acceptable; > 0.2: Potentially misleading
cat("NMDS Stress Value:", round(nmds_res$stress, 4), "\n")

# 3. Create a dataframe for plotting
nmds_data <- data.frame(
  NMDS1   = nmds_res$points[, 1],
  NMDS2   = nmds_res$points[, 2],
  Cluster = factor(cluster_sub)
)
my_palette <- c("#93B5C6", "#DDEDAA", "#D7816A", "#F0CF65", "#BD4F6C")
# 4. Visualization with ggplot2
nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_jitter(width = 0.3, height = 0.3, alpha = 0.8, size = 2) +
  scale_color_manual(values = my_palette) +
  theme_minimal() +
  labs(
    title = "Non-metric MDS: Respondents clusters",
    x = "Dimension 1",
    y = " Dimension 2"
  ) +
  theme(plot.title = element_text(face = "bold"))

nmds_plot
