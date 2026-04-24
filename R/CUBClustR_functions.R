# library(Rcpp)
# library(mclust)
# library(ggplot2)
# library(reshape2)
# library(doParallel)
# library(foreach)

# ==============================================================================
# INTERNAL PATH DETECTION (HIDDEN)
# ==============================================================================
#
# .get_script_path <- function() {
#   cmdArgs <- commandArgs(trailingOnly = FALSE)
#   needle <- "--file="
#   match <- grep(needle, cmdArgs)
#
#   if (length(match) > 0) {
#     return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
#   } else {
#     # Detect path if sourced in RStudio or standard R GUI
#     filename <- NULL
#     try({ filename <- sys.frames()[[1]]$filename }, silent = TRUE)
#
#     if (!is.null(filename)) {
#       return(dirname(normalizePath(filename)))
#     }
#     return(getwd())
#   }
# }

# Use the hidden function to locate and compile the C++ engine
# .script_dir <- .get_script_path()
# .cpp_path   <- file.path(.script_dir, "engine.cpp")
#
# if (file.exists(.cpp_path)) {
#   sourceCpp(.cpp_path)
# } else {
#   warning("C++ engine (engine.cpp) not found in the script directory.")
# }

# ==============================================================================
# 1. FIT MODEL
# ==============================================================================

#' Fit Multivariate Latent Class CUB Model
#'
#' @description Estimates the parameters of a Multivariate Latent Class CUB model using an EM algorithm accelerated by C++.
#'
#' @param data A numeric matrix or data.frame containing the observed ordinal ratings.
#' @param m A numeric vector specifying the maximum number of categories for each item.
#' @param K The number of latent classes (clusters) to estimate.
#' @param tol Tolerance for the log-likelihood convergence criterion (default: 1e-5).
#' @param EM_iter Number of random initializations for the EM algorithm to avoid local maxima (default: 10).
#' @param max_iter Maximum number of iterations allowed for a single EM run (default: 1500).
#'
#' @return A list containing the estimated model parameters (`pi`, `xi`, `omega`), the final log-likelihood (`LL`), the assigned cluster for each subject (`class`), and the number of iterations until convergence (`iterations`).
#' @examples
#' \dontrun{
#' # Generate a mock dataset (100 respondents, 4 items, 5-point scale)
#' set.seed(123)
#' data_mock <- matrix(sample(1:5, 100 * 4, replace = TRUE), ncol = 4)
#' m_vec <- rep(5, 4)
#' 
#' # Fit the model with 2 clusters
#' # (Using very few iterations for demonstration purposes only)
#' model <- fit_MLCCUB(data = data_mock, m = m_vec, K = 2, EM_iter = 2, max_iter = 50)
#' 
#' # Display the final log-likelihood and cluster mixing proportions
#' print(model$LL)
#' print(model$omega)
#' }
#' @export
fit_MLCCUB <- function(data, m, K, tol = 1e-5, EM_iter = 10, max_iter = 1500) {

  R <- as.matrix(data)
  storage.mode(R) <- "double"
  m <- as.numeric(m)

  n <- nrow(R)
  J <- ncol(R)

  uniform.j <- matrix(NA, nrow = n, ncol = J)
  for(j in 1:J) {
    uniform.j[, j] <- rep(1 / m[j], n)
  }
  storage.mode(uniform.j) <- "double"

  best_model <- NULL
  best_ll <- -Inf

  for(em in 1:EM_iter) {

    pi_mat_curr <- matrix(runif(K * J, min = 0.1, max = 0.9), nrow = K, ncol = J)
    xi_mat_curr <- matrix(runif(K * J, min = 0.1, max = 0.9), nrow = K, ncol = J)
    omega_curr <- rep(1 / K, K)

    ll_prev <- -Inf
    iter <- 0

    success <- tryCatch({

      while(iter < max_iter) {

        # E-step (C++)
        LL_out <- .LogLik_cpp(as.numeric(omega_curr), pi_mat_curr, R, m, xi_mat_curr, uniform.j)
        tau_mat <- sweep(LL_out$multiCUB_w_k, 1, LL_out$MLCCUBk + 1e-15, "/")

        # Reshape eta flat vector from C++
        eta_tens_flat <- .eta_ijk_cpp(pi_mat_curr, R, m, xi_mat_curr, uniform.j, K, J, n)
        dim(eta_tens_flat) <- c(K, J, n)
        eta_tens <- eta_tens_flat

        # M-step
        omega_curr <- colMeans(tau_mat)

        for(k in 1:K) {
          sum_tau_k <- sum(tau_mat[, k]) + 1e-15

          for(j in 1:J) {
            pi_new <- sum(tau_mat[, k] * eta_tens[k, j, ]) / sum_tau_k
            pi_mat_curr[k, j] <- min(max(pi_new, 1e-4), 1 - 1e-4)

            num_xi <- sum(tau_mat[, k] * eta_tens[k, j, ] * (m[j] - R[, j]))
            den_xi <- sum(tau_mat[, k] * eta_tens[k, j, ] * (m[j] - 1))
            xi_new <- num_xi / (den_xi + 1e-15)
            xi_mat_curr[k, j] <- min(max(xi_new, 1e-4), 1 - 1e-4)
          }
        }

        new_ll <- LL_out$LL
        if(is.finite(new_ll) && is.finite(ll_prev)) {
          if(abs(new_ll - ll_prev) < tol) break
        }
        ll_prev <- new_ll
        iter <- iter + 1
      }
      TRUE
    }, error = function(e) {
      message("    [Debug] EM iter ", em, " failed: ", e$message)
      FALSE
    })

    if(success && is.finite(ll_prev) && ll_prev > best_ll) {
      best_ll <- ll_prev
      final_LL_out <- .LogLik_cpp(as.numeric(omega_curr), pi_mat_curr, R, m, xi_mat_curr, uniform.j)
      final_tau <- sweep(final_LL_out$multiCUB_w_k, 1, final_LL_out$MLCCUBk + 1e-15, "/")

      best_model <- list(
        pi = pi_mat_curr,
        xi = xi_mat_curr,
        omega = omega_curr,
        LL = best_ll,
        class = max.col(final_tau),
        iterations = iter
      )
    }
  }

  if(is.null(best_model)) {
    stop("The algorithm failed to converge. Check your data or increase EM_iter/tolerance.")
  }

  return(best_model)
}

# ==============================================================================
# 2. COMPARE MODELS
# ==============================================================================

#' Compare MLCCUB Models for Different Values of K
#'
#' @description Fits multiple MLCCUB models for a sequence of latent classes (K) and compares them using the Bayesian Information Criterion (BIC).
#'
#' @param data A numeric matrix or data.frame containing the observed ordinal ratings.
#' @param m A numeric vector specifying the maximum number of categories for each item.
#' @param K_vector A numeric vector containing the different values of K to be tested (e.g., c(2, 3, 4)).
#' @param ... Additional arguments passed to the `fit_MLCCUB` function.
#'
#' @return A list containing all the fitted models (`all_models`), a summary data.frame with LogLikelihood and BIC (`summary`), the best K value (`best_k`), and the best model object (`best_model`).
#' @examples
#' \dontrun{
#' set.seed(123)
#' data_mock <- matrix(sample(1:7, 100 * 3, replace = TRUE), ncol = 3)
#' m_vec <- rep(7, 3)
#' 
#' # Quick comparison among models with K = 1, 2, and 3
#' comparison <- compare_MLCCUB_models(data_mock, m_vec, K_vector = 1:3, EM_iter = 2)
#' 
#' # Print the summary with LogLik and BIC to choose the optimal model
#' print(comparison$summary)
#' }
#' @export
compare_MLCCUB_models <- function(data, m, K_vector, ...) {

  n_subjects <- nrow(data)
  J_items <- ncol(data)

  all_models <- list()
  summary_stats <- data.frame(
    K = K_vector,
    LogLik = NA,
    n_params = NA,
    BIC = NA
  )

  for (i in seq_along(K_vector)) {
    k_val <- K_vector[i]
    message("Analyzing K = ", k_val, "...")

    model_fit <- tryCatch({
      fit_MLCCUB(data = data, m = m, K = k_val, ...)
    }, error = function(e) {
      message("  -> Error for K = ", k_val, ": ", e$message)
      return(NULL)
    })

    if (!is.null(model_fit)) {
      n_p <- (k_val * J_items * 2) + (k_val - 1)
      bic_val <- -2 * model_fit$LL + n_p * log(n_subjects)

      summary_stats$LogLik[i] <- model_fit$LL
      summary_stats$n_params[i] <- n_p
      summary_stats$BIC[i] <- bic_val

      model_fit$K <- k_val
      model_fit$BIC <- bic_val

      all_models[[paste0("K", k_val)]] <- model_fit
    }
  }

  best_idx <- which.min(summary_stats$BIC)
  best_k <- summary_stats$K[best_idx]

  results <- list(
    all_models = all_models,
    summary = summary_stats,
    best_k = best_k,
    best_model = all_models[[paste0("K", best_k)]]
  )

  message("\nBest model selected according to BIC: K = ", best_k)
  return(results)
}

# ==============================================================================
# 3. PLOT PARAMETERS
# ==============================================================================

#' Plot MLCCUB Parameters in the Parameter Space
#'
#' @description Generates a 2D plot showing the Feeling and Uncertainty parameters for each item, grouped by latent cluster.
#'
#' @param model_obj The model object returned by `fit_MLCCUB`.
#' @param cluster_order Optional numeric vector to manually reorder the clusters (e.g., c(1, 3, 2)).
#' @param custom_colors Optional character vector specifying custom HEX colors for the clusters.
#'
#' @return A `ggplot2` object displaying the parameter space.
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_point geom_text scale_color_manual scale_fill_manual theme_classic ylab xlab theme element_text unit element_blank
#' @examples
#' \dontrun{
#' set.seed(123)
#' data_mock <- matrix(sample(1:5, 100 * 4, replace = TRUE), ncol = 4)
#' model <- fit_MLCCUB(data_mock, m = rep(5, 4), K = 2, EM_iter = 1, max_iter = 20)
#' 
#' # Create the parameter space plot (Feeling vs. Uncertainty)
#' p_space <- plot_MLCCUB_clusters(model)
#' print(p_space)
#' }
#' @export
plot_MLCCUB_clusters <- function(model_obj, cluster_order = NULL, custom_colors = NULL) {

  pi_mat <- model_obj$pi
  xi_mat <- model_obj$xi
  omega_vec <- model_obj$omega

  K <- length(omega_vec)
  J <- ncol(pi_mat)

  if (!is.null(cluster_order)) {
    if (length(cluster_order) != K) {
      stop("The length of the 'cluster_order' vector must be equal to K.")
    }
    pi_mat <- pi_mat[cluster_order, , drop = FALSE]
    xi_mat <- xi_mat[cluster_order, , drop = FALSE]
    omega_vec <- omega_vec[cluster_order]
  }

  x_df <- as.data.frame(round(1 - pi_mat, 2))
  y_df <- as.data.frame(round(1 - xi_mat, 2))

  colnames(x_df) <- as.character(1:J)
  colnames(y_df) <- as.character(1:J)

  x_df$cluster <- 1:K
  y_df$cluster <- 1:K

  x_long <- melt(x_df, id.vars = "cluster", variable.name = "variable", value.name = "x")
  y_long <- melt(y_df, id.vars = "cluster", variable.name = "variable", value.name = "y")

  x_long$variable_number <- as.numeric(as.character(x_long$variable))
  y_long$variable_number <- as.numeric(as.character(y_long$variable))

  plot_data <- merge(x_long, y_long, by = c("cluster", "variable", "variable_number"))

  omega_round <- round(omega_vec, 2)
  legend_labels <- lapply(1:K, function(i) {
    bquote(atop("Cluster"~.(i), ~omega[.(i)]~"="~.(omega_round[i])))
  })
  legend_labels <- lapply(legend_labels, as.expression)

  if (!is.null(custom_colors)) {
    base_colors <- custom_colors
  } else {
    base_colors <- c("#93B5C6", "#DDEDAA", "#D7816A", "#F0CF65", "#BD4F6C", "#7C6A92", "#ECA784")
  }

  if (K <= length(base_colors)) {
    plot_colors <- base_colors[1:K]
  } else {
    plot_colors <- colorRampPalette(base_colors)(K)
  }

  p <- ggplot(plot_data, aes(x = x, y = y, color = factor(cluster), fill = factor(cluster))) +
    geom_point(shape = 21, size = 5) +
    geom_text(aes(label = variable_number), color = "black", size = 4, vjust = 0.5, hjust = 0.5) +
    scale_color_manual(name = "Cluster", values = plot_colors, labels = legend_labels) +
    scale_fill_manual(name = "Cluster", values = plot_colors, labels = legend_labels) +
    theme_classic() +
    ylab(expression("Feeling (1 - " * xi[jk] * ")")) +
    xlab(expression("Uncertainty (1 - " * pi[jk] * ")")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))

  return(p)
}

# ==============================================================================
# 4. PLOT FIT (OBSERVED VS THEORETICAL)
# ==============================================================================

#' Minimal Plot of CUB Distributions (Observed vs Theoretical)
#'
#' @description Plots the observed relative frequencies of the responses against the theoretical probability distributions estimated by the MLCCUB model.
#'
#' @param model_obj The model object returned by `fit_MLCCUB`.
#' @param data The original dataset containing the observed responses.
#' @param m A numeric vector specifying the maximum number of categories for each item.
#' @param cluster_order Optional numeric vector to manually reorder the clusters.
#' @param line_color Color of the dashed theoretical probability line (default: "#008B8B").
#'
#' @return A `ggplot2` object showing the fit comparison.
#' @importFrom ggplot2 ggplot aes geom_segment geom_line geom_point facet_grid scale_x_continuous scale_color_manual theme_bw labs theme element_blank element_text unit
#' @importFrom stats dbinom
#' @examples
#' \dontrun{
#' set.seed(123)
#' data_mock <- matrix(sample(1:5, 100 * 4, replace = TRUE), ncol = 4)
#' m_vec <- rep(5, 4)
#' model <- fit_MLCCUB(data_mock, m = m_vec, K = 2, EM_iter = 1, max_iter = 20)
#' 
#' # Goodness-of-fit plot: theoretical probabilities vs observed frequencies
#' p_fit <- plot_CUB_fit(model, data = data_mock, m = m_vec)
#' print(p_fit)
#' }
#' @export
plot_CUB_fit <- function(model_obj, data, m, cluster_order = NULL, line_color = "darkblue") {

  pi_mat <- model_obj$pi
  xi_mat <- model_obj$xi
  omega_vec <- model_obj$omega
  assignments <- model_obj$class

  K <- length(omega_vec)
  J <- ncol(pi_mat)

  # Assicuriamoci che data sia una matrice numerica
  R_data <- as.matrix(data)

  # --- DEBUG CHECK ---
  if (max(assignments) > nrow(R_data)) {
    stop(paste("Error: The model assignments refer to index", max(assignments),
               "but the data provided only has", nrow(R_data), "rows."))
  }
  # --------------------

  if (length(m) != J) stop("The length of vector 'm' must equal the number of items.")

  plot_data <- data.frame()

  for (k in 1:K) {
    idx_k <- which(assignments == k)

    # Se il cluster è vuoto, lo saltiamo
    if (length(idx_k) == 0) next

    for (j in 1:J) {
      m_j <- m[j]
      r_seq <- 1:m_j

      pi_kj <- pi_mat[k, j]
      xi_kj <- xi_mat[k, j]

      # Distribuzione teorica CUB
      prob_theo <- pi_kj * dbinom(r_seq - 1, size = m_j - 1, prob = 1 - xi_kj) + (1 - pi_kj) * (1 / m_j)

      # Estrazione sicura delle risposte osservate per il cluster k e item j
      risposte_item <- R_data[idx_k, j]
      risposte_item <- risposte_item[!is.na(risposte_item)]

      if (length(risposte_item) > 0) {
        counts <- table(factor(risposte_item, levels = r_seq))
        prob_obs <- as.numeric(counts) / length(risposte_item)
      } else {
        prob_obs <- rep(0, m_j)
      }

      plot_data <- rbind(plot_data, data.frame(
        Cluster_Orig = k,
        Item = paste0("Item ", j),
        Category = r_seq,
        Prob_Theo = prob_theo,
        Prob_Obs = prob_obs
      ))
    }
  }

  # Gestione ordinamento cluster
  if (!is.null(cluster_order)) {
    plot_data$Cluster_Num <- match(plot_data$Cluster_Orig, cluster_order)
  } else {
    plot_data$Cluster_Num <- plot_data$Cluster_Orig
  }

  plot_data$Cluster <- factor(paste0("Cluster ", plot_data$Cluster_Num),
                              levels = paste0("Cluster ", 1:K))
  plot_data$Item <- factor(plot_data$Item, levels = paste0("Item ", 1:J))

  # Creazione del grafico con ggplot2
  p <- ggplot(plot_data, aes(x = Category)) +
    geom_segment(aes(xend = Category, y = 0, yend = Prob_Obs, color = "Observed Frequency"), linewidth = 1) +
    geom_line(aes(y = Prob_Theo, color = "Theoretical Probability"), linetype = "dashed", linewidth = 0.8) +
    geom_point(aes(y = Prob_Theo, color = "Theoretical Probability"), size = 2.5) +
    facet_grid(Cluster ~ Item) +
    scale_x_continuous(breaks = 1:max(m)) +
    scale_color_manual(
      name = "",
      values = c("Observed Frequency" = "black", "Theoretical Probability" = line_color)
    ) +
    theme_bw() +
    labs(x = "Response Category", y = "Probability / Relative Frequency") +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )

  return(p)
}
# ==============================================================================
# 5. BOOTSTRAP IDENTIFIABILITY
# ==============================================================================

#' Model Identifiability Check via Bootstrap
#'
#' @description Performs a non-parametric bootstrap to assess the stability and identifiability of the latent clusters using the Adjusted Rand Index (ARI).
#'
#' @param data The original dataset (matrix or data.frame).
#' @param m A numeric vector specifying the maximum number of categories for each item.
#' @param K The number of clusters to test.
#' @param n_boot Number of bootstrap samples to extract (default: 100).
#' @param n_cores Number of CPU cores for parallel computation (default: 1 for sequential).
#' @param plot_hist Logical. If TRUE, generates and returns the ARI distribution histogram.
#' @param seed Random seed for reproducibility.
#' @param cpp_file Path to the C++ engine file to be loaded by workers (default: "engine.cpp").
#' @param ... Additional arguments passed to `fit_MLCCUB` (e.g., `EM_iter`, `tol`).
#'
#' @return A list containing the pairwise ARI matrix (`ARI_matrix`), the flattened ARI values (`ARI_values`), the ggplot object (`plot`), and the raw cluster assignments (`raw_clusters`).
#' @importFrom parallel makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom ggplot2 ggplot geom_histogram labs theme_classic theme element_text xlim aes
#' @examples
#' \dontrun{
#' set.seed(123)
#' data_mock <- matrix(sample(1:5, 100 * 4, replace = TRUE), ncol = 4)
#' 
#' # Fast bootstrap with only 5 iterations and 1 core (for testing purposes)
#' ident_res <- check_identifiability(data = data_mock, m = rep(5, 4), K = 2, 
#'                                    n_boot = 5, n_cores = 1, EM_iter = 1)
#' 
#' # Display the histogram of the Adjusted Rand Index (ARI) distribution
#' print(ident_res$plot)
#' }
#' @export
check_identifiability <- function(data, m, K, n_boot = 100, n_cores = 1,
                                  plot_hist = TRUE, seed = 250535,
                                  cpp_file = .cpp_path, ...) {

  set.seed(seed)
  n <- nrow(data)

  if (n_cores > 1) {
    message("Configuring parallel cluster with ", n_cores, " cores...")
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    clusterExport(cl, varlist = c("fit_MLCCUB", "cpp_file"), envir = environment())
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_file)
    })
  } else {
    registerDoSEQ()
  }

  message("Starting Bootstrap (", n_boot, " iterations) for K = ", K, "...")

  boot_results <- foreach(b = 1:n_boot, .packages = c("Rcpp")) %dopar% {

    boot_idx <- sample(1:n, replace = TRUE)
    boot_data <- data[boot_idx, ]

    mod <- fit_MLCCUB(data = boot_data, m = m, K = K, ...)

    list(index = boot_idx, class = mod$class)
  }

  if (n_cores > 1) {
    stopCluster(cl)
  }

  message("Computing pairwise ARI matrix...")
  num_elements <- length(boot_results)
  ARI_matrix <- matrix(NA, nrow = num_elements, ncol = num_elements)

  for (i in 1:num_elements) {
    index1 <- boot_results[[i]]$index
    for (j in i:num_elements) {
      index2 <- boot_results[[j]]$index

      intersection <- intersect(index1, index2)

      if(length(intersection) > 0) {
        match1 <- match(intersection, index1)
        match2 <- match(intersection, index2)

        class1 <- boot_results[[i]]$class[match1]
        class2 <- boot_results[[j]]$class[match2]

        ARI_val <- mclust::adjustedRandIndex(class1, class2)
      } else {
        ARI_val <- NA
      }

      ARI_matrix[i, j] <- ARI_val
      ARI_matrix[j, i] <- ARI_val
    }
  }

  diag(ARI_matrix) <- NA
  ARI_vector <- ARI_matrix[upper.tri(ARI_matrix)]
  ARI_vector <- ARI_vector[!is.na(ARI_vector)]

  p <- NULL
  if (plot_hist) {
    ARI_df <- data.frame(ARI = ARI_vector)
    p <- ggplot(ARI_df, aes(x = ARI)) +
      geom_histogram(binwidth = 0.03, color = "black", fill = "lightgrey", alpha = 0.8) +
      labs(x = "Adjusted Rand Index (ARI)",
           y = "Frequency",
           title = paste("Bootstrap ARI Distribution (K =", K, ")")) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
      xlim(-1, 1)
  }

  message("Procedure completed!")

  return(list(
    ARI_matrix = ARI_matrix,
    ARI_values = ARI_vector,
    plot = p,
    raw_clusters = boot_results
  ))
}

# ==============================================================================
# 6. CALCULATE WEIGHTED AVERAGE DISSIMILARITY (WAD) AND DIAGNOSTICS
# ==============================================================================

#' Calculate Weighted Average Dissimilarity (WAD) and Detailed Dissimilarities
#'
#' @description Computes the WAD index and detailed dissimilarities to assess the
#' goodness of fit of the MLC-CUB model. Returns values in `[0, 1]`, where 0
#' indicates a perfect fit between observed data and theoretical probabilities.
#'
#' @param model_obj The model object returned by `fit_MLCCUB`.
#' @param data The original dataset containing the observed ordinal responses.
#' @param m A numeric vector specifying the maximum number of categories for each item.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{WAD}: Global Weighted Average Dissimilarity index.
#'   \item \code{Cluster_Dissimilarity}: Unweighted average dissimilarity for each cluster.
#'   \item \code{Item_Dissimilarity}: Weighted average dissimilarity for each item.
#'   \item \code{Matrix_Dissimilarity}: A K x J matrix of normalized dissimilarities.
#' }
#' @importFrom stats dbinom
#' @examples
#' \dontrun{
#' set.seed(123)
#' data_mock <- matrix(sample(1:5, 100 * 4, replace = TRUE), ncol = 4)
#' m_vec <- rep(5, 4)
#' model <- fit_MLCCUB(data_mock, m = m_vec, K = 2, EM_iter = 1, max_iter = 20)
#' 
#' # Calculate the Weighted Average Dissimilarity (WAD) indices
#' wad_metrics <- calculate_WAD(model_obj = model, data = data_mock, m = m_vec)
#' 
#' # Print the global WAD index (values closer to 0 indicate a better fit)
#' print(wad_metrics$WAD)
#' }
#' @export
calculate_WAD <- function(model_obj, data, m) {

  pi_mat <- model_obj$pi
  xi_mat <- model_obj$xi
  omega_vec <- model_obj$omega
  assignments <- model_obj$class

  K <- length(omega_vec)
  J <- ncol(pi_mat)

  # Robust matrix conversion to avoid "subscript out of bounds"
  R_data <- as.matrix(data)
  storage.mode(R_data) <- "numeric"

  # Integrity check: Ensure data rows match model assignment length
  if (max(assignments) > nrow(R_data)) {
    stop(paste("Error: The model assignments refer to index", max(assignments),
               "but the data provided only has", nrow(R_data), "rows."))
  }

  if (length(m) != J) stop("The length of vector 'm' must equal the number of items.")

  # Initialize matrix to store normalized [0,1] dissimilarities for each cluster-item pair
  diss_matrix <- matrix(0, nrow = K, ncol = J)
  rownames(diss_matrix) <- paste0("Cluster_", 1:K)
  colnames(diss_matrix) <- paste0("Item_", 1:J)

  for (j in 1:J) {
    m_j <- m[j]
    r_seq <- 1:m_j

    for (k in 1:K) {
      idx_k <- which(assignments == k)

      # 1. Calculate Observed Frequencies for Cluster k, Item j
      if (length(idx_k) > 0) {
        item_responses <- R_data[idx_k, j]
        item_responses <- item_responses[!is.na(item_responses)]

        if (length(item_responses) > 0) {
          counts <- table(factor(item_responses, levels = r_seq))
          f_obs <- as.numeric(counts) / length(item_responses)
        } else {
          f_obs <- rep(0, m_j)
        }
      } else {
        f_obs <- rep(0, m_j)
      }

      # 2. Calculate Theoretical Expected Probabilities
      pi_kj <- pi_mat[k, j]
      xi_kj <- xi_mat[k, j]
      p_exp <- pi_kj * dbinom(r_seq - 1, size = m_j - 1, prob = 1 - xi_kj) + (1 - pi_kj) * (1 / m_j)

      # 3. Compute Normalized Dissimilarity (Total Variation Distance)
      # Divided by 2 to bound the value within [0, 1]
      diss_jk <- sum(abs(f_obs - p_exp)) / 2
      diss_matrix[k, j] <- diss_jk
    }
  }

  # A. Cluster-specific Dissimilarity (Unweighted average across items)
  cluster_diss <- rowMeans(diss_matrix)

  # B. Item-specific Dissimilarity (Average weighted by cluster sizes omega_k)
  item_diss <- colSums(diss_matrix * omega_vec)

  # C. Global WAD (Weighted mean of all dissimilarities)
  wad_value <- mean(item_diss)

  return(list(
    WAD = wad_value,
    Cluster_Dissimilarity = cluster_diss,
    Item_Dissimilarity = item_diss,
    Matrix_Dissimilarity = diss_matrix
  ))
}
