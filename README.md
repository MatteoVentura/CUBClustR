# CUBClustR: Clustering Rating Data within the CUB Framework

**CUBClustR** is an R package that implements the **MLC-CUB (Multivariate Latent Class CUB)** model. It allows researchers to segment ordinal data (ratings) by accounting for both the *feeling process* (reasoned choice) and the *uncertainty* of the respondent, grouping subjects into $K$ latent clusters. 

The estimation relies on the Expectation-Maximization (EM) algorithm, deeply optimized via a **C++ backend (Rcpp)** to ensure high computational performance.

## Installation

You can install the development version of `CUBClustR` directly from GitHub. Ensure you have the `vegan` package installed first:

```R
# Install required dependencies
install.packages("vegan")

# Install CUBClustR
devtools::install_github("MatteoVentura/CUBClustR", build_vignettes = FALSE)
```
Quick Start (University Satisfaction Dataset)
Here is a quick example of how to use CUBClustR to analyze student satisfaction across 5 items (7-point scale) using the univer dataset from the CUB package.

## 1. Setup and Data Preparation
```R
library(CUBClustR)
library(CUB)
library(ggplot2)

# Load data and select 5 relevant ordinal items
data(univer)
R <- apply(univer[, 8:12], 2, as.integer)

# Maximum categories for each item (m = 7)
m_vector <- rep(7, ncol(R))
```

## 2. Model Selection (BIC)
Test from 1 to 7 latent classes to find the optimal balance between model complexity and fit:

```R
model_selection <- compare_MLCCUB_models(
  data     = R, 
  m        = m_vector, 
  K_vector = 1:7, 
  EM_iter  = 10, 
  max_iter = 1000
)
```

### Plotting the BIC trend
```R
ggplot(model_selection$summary, aes(x = K, y = BIC)) +
  geom_line(color = "#2C3E50", linewidth = 1.2) +
  geom_point(color = "#BD4F6C", size = 3) +
  theme_bw() +
  labs(title = "Model Selection: BIC Trend (University Data)", x = "Number of Latent Classes (K)")
```

## 3. Stability & Identifiability Check
Evaluate the identifiability of a 5-cluster solution using parallelized bootstrap:

```R
identifiability <- check_identifiability(
  data      = R, 
  m         = m_vector, 
  K         = 5, 
  n_boot    = 100, 
  n_cores   = 2  # Adjust based on your CPU
)
```
Display the ARI distribution: values near 1 indicate highly stable clusters
```R
print(identifiability$plot)
```

## 4. Visualization & Goodness of Fit
Visualize the parameter space and calculate the absolute fit quality (WAD index):

```R
model_k5 <- model_selection$all_models$K5
```

### Parameter Space Plot (Feeling vs. Uncertainty)
```R
plot_MLCCUB_clusters(model_k5) + 
  labs(title = "University Clusters: Parameter Space")
```

### Goodness of Fit metrics
```R
fit_metrics <- calculate_WAD(model_obj = model_k5, data = R, m = m_vector)
print(fit_metrics$WAD) # Global index (closer to 0 is better)
```

📚 Documentation
For a complete overview of the mathematical framework, diagnostic indices, and Non-metric Multidimensional Scaling (NMDS) examples, please check the following resources located in the docs/ folder:

* **[Official Technical Manual (PDF)](docs/CUBClustR_manual.pdf)**
* **[Full Example Script (R)](docs/univ_script.R)**
