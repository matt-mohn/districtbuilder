# Package Load ------------------------------------------------------------

library(sf)
library(sfdep)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

library(sf)
library(dplyr)
library(cluster)

initialize_districts_from_shapefile <- function(shapefile_path,
                                                d,
                                                clustering_method = "spatial",
                                                seed = 123,
                                                geoid_column = "geoid",
                                                verbose = TRUE) {
  set.seed(seed)
  if (verbose) cat("📁 Loading shapefile:", shapefile_path, "\n")
  
  # Load shapefile
  if (!file.exists(shapefile_path)) {
    stop("Shapefile not found: ", shapefile_path)
  }
  
  nc_data <- read_sf(shapefile_path)
  if (verbose) cat("✅ Loaded", nrow(nc_data), "features\n")
  
  # Parse attributes as numeric ---------------------------------------------
  if (verbose) cat("🔢 Converting attributes to numeric (preserving", geoid_column, ")\n")
  
  # Get column names
  all_cols <- names(nc_data)
  geom_col <- attr(nc_data, "sf_column") # geometry column name
  
  # Identify columns to convert (exclude geoid and geometry)
  cols_to_convert <- all_cols[!all_cols %in% c(geoid_column, geom_col)]
  
  # Convert to numeric where possible
  for (col in cols_to_convert) {
    original_class <- class(nc_data[[col]])[1]
    # Try to convert to numeric
    numeric_version <- suppressWarnings(as.numeric(as.character(nc_data[[col]])))
    # Check if conversion was successful (no NAs introduced)
    if (sum(is.na(numeric_version)) <= sum(is.na(nc_data[[col]]))) {
      nc_data[[col]] <- numeric_version
      if (verbose) cat("  ✓", col, ":", original_class, "->", "numeric\n")
    } else {
      if (verbose) cat("  ⚠", col, ":", original_class, "-> keeping original (non-numeric)\n")
    }
  }
  
  # Spatial clustering ------------------------------------------------------
  if (verbose) cat("🗺️ Creating", d, "districts using", clustering_method, "clustering\n")
  
  if (clustering_method == "spatial") {
    # Use spatial coordinates of centroids
    centroids <- st_centroid(nc_data)
    coords <- st_coordinates(centroids)
    # Standardize coordinates
    coords_scaled <- scale(coords)
    # K-means clustering on coordinates
    clusters <- kmeans(coords_scaled, centers = d, nstart = 25, iter.max = 100)
    nc_data$district <- as.integer(clusters$cluster)
  } else if (clustering_method == "kmeans_coords") {
    # Similar to spatial but without scaling
    centroids <- st_centroid(nc_data)
    coords <- st_coordinates(centroids)
    clusters <- kmeans(coords, centers = d, nstart = 25, iter.max = 100)
    nc_data$district <- as.integer(clusters$cluster)
  } else if (clustering_method == "kmeans_attributes") {
    # Use numeric attributes for clustering (excluding spatial info)
    numeric_cols <- sapply(nc_data, is.numeric)
    numeric_cols[geom_col] <- FALSE # Exclude geometry
    if (sum(numeric_cols) == 0) {
      stop("No numeric columns found for attribute-based clustering")
    }
    # Extract numeric data
    numeric_data <- st_drop_geometry(nc_data[, numeric_cols, drop = FALSE])
    # Remove any rows with NAs
    complete_rows <- complete.cases(numeric_data)
    if (sum(complete_rows) < nrow(numeric_data)) {
      if (verbose) cat("  ⚠ Removing", sum(!complete_rows), "rows with missing data\n")
      nc_data <- nc_data[complete_rows, ]
      numeric_data <- numeric_data[complete_rows, ]
    }
    # Standardize data
    numeric_data_scaled <- scale(numeric_data)
    # K-means clustering
    clusters <- kmeans(numeric_data_scaled, centers = d, nstart = 25, iter.max = 100)
    nc_data$district <- as.integer(clusters$cluster)
  } else {
    stop("Unknown clustering method: ", clustering_method)
  }
  
  # Pre-calculate geometric properties for fast evaluation ------------------
  if (verbose) cat("🚀 Pre-calculating geometric properties for lightning-fast evaluation...\n")
  
  # Calculate centroids for each precinct (needed for bounding box compactness)
  centroids <- st_centroid(nc_data)
  coords <- st_coordinates(centroids)
  nc_data$centroid_x <- coords[, 1]
  nc_data$centroid_y <- coords[, 2]
  
  if (verbose) cat("✅ Pre-calculation complete!\n")
  
  # Print summary statistics if verbose -------------------------------------
  if (verbose) {
    cat("✅ District initialization complete!\n\n")
    cat("📊 Summary:\n")
    cat("  Total features:", nrow(nc_data), "\n")
    cat("  Districts created:", d, "\n")
    cat("  Clustering method:", clustering_method, "\n")
    
    # District sizes
    district_counts <- table(nc_data$district)
    cat("  District sizes:\n")
    for (i in 1:d) {
      cat("    District", i, ":", district_counts[i], "precincts\n")
    }
    
    # Column summary
    cat("  Columns in dataset:", ncol(nc_data), "\n")
    cat("    Numeric columns:", sum(sapply(nc_data, is.numeric)), "\n")
    cat("    Character columns:", sum(sapply(nc_data, is.character)), "\n")
    cat("    Geometry column: 1\n")
  }
  
  return(nc_data)
}

# Function: Calculate contiguity score using connected components
calculate_contiguity_score <- function(nc_data, neighbours_queen) {
  districts <- unique(nc_data$district)
  n_districts <- length(districts)
  total_components <- 0
  
  for (d in districts) {
    district_precincts <- which(nc_data$district == d)
    
    if (length(district_precincts) == 0) {
      next
    }
    
    if (length(district_precincts) == 1) {
      total_components <- total_components + 1
      next
    }
    
    # Create adjacency list for this district only
    district_adj <- vector("list", length(district_precincts))
    names(district_adj) <- as.character(district_precincts)
    
    for (i in seq_along(district_precincts)) {
      precinct_idx <- district_precincts[i]
      # Get neighbors that are also in this district
      neighbors_in_district <- intersect(neighbours_queen[[precinct_idx]], district_precincts)
      # Convert to local indices within this district
      local_neighbors <- which(district_precincts %in% neighbors_in_district)
      district_adj[[i]] <- local_neighbors
    }
    
    # Count connected components using DFS
    visited <- rep(FALSE, length(district_precincts))
    components <- 0
    
    dfs <- function(node) {
      visited[node] <<- TRUE
      for (neighbor in district_adj[[node]]) {
        if (!visited[neighbor]) {
          dfs(neighbor)
        }
      }
    }
    
    for (i in seq_along(district_precincts)) {
      if (!visited[i]) {
        components <- components + 1
        dfs(i)
      }
    }
    
    total_components <- total_components + components
  }
  
  # Calculate contiguity score as percentage  
  # Perfect: total_components = n_districts (1.0 = 100% contiguity)
  # Realistic worst case: total_components = n_districts * 10 (districts broken into ~10 pieces each)
  worst_case <- n_districts * 10  # Realistic worst case
  perfect_case <- n_districts
  
  if (total_components <= perfect_case) {
    contiguity_score <- 1.0  # Perfect contiguity
  } else if (total_components >= worst_case) {
    contiguity_score <- 0.0  # At or beyond worst case
  } else {
    # Linear scale between perfect and worst case
    contiguity_score <- 1.0 - (total_components - perfect_case) / (worst_case - perfect_case)
  }
  
  return(list(
    score = contiguity_score,
    total_components = total_components,
    target_components = n_districts
  ))
}

# Function: Calculate border precincts (uses rook contiguity)
calculate_border <- function(nc_data, neighbours_rook) {
  border <- rep(0, nrow(nc_data))
  neighboring_districts <- rep("", nrow(nc_data))
  
  for (i in seq_along(neighbours_rook)) {
    current_district <- nc_data$district[i]
    neighbor_districts <- unique(nc_data$district[neighbours_rook[[i]]])
    other_districts <- neighbor_districts[neighbor_districts != current_district]
    
    if (length(other_districts) > 0) {
      border[i] <- 1
      neighboring_districts[i] <- paste(sort(other_districts), collapse = ",")
    }
  }
  
  return(list(border = border, neighboring_districts = neighboring_districts))
}

# Function: Update analysis for current district assignment
update_analysis <- function(nc_data, neighbours_queen, neighbours_rook) {
  border_results <- calculate_border(nc_data, neighbours_rook)
  nc_data$border <- border_results$border
  nc_data$neighboring_districts <- border_results$neighboring_districts
  return(nc_data)
}

# Function: Execute multiple swaps (1 to s swaps per edit)
execute_edit <- function(nc_data, neighbours_queen, neighbours_rook, s) {
  current_data <- nc_data
  
  # Update analysis before swaps
  current_data <- update_analysis(current_data, neighbours_queen, neighbours_rook)
  
  # Randomly choose number of swaps between 1 and s (biased toward lower values)
  weights <- rev(2^(0:(s - 1))) # Creates weights favoring fewer swaps
  actual_s <- sample(1:s, 1, prob = weights)
  
  swaps_made <- 0
  
  for (i in 1:actual_s) {
    # Get eligible border precincts
    swapworthy <- which(current_data$border == 1)
    
    if (length(swapworthy) == 0) {
      break
    }
    
    # Safe sampling
    if (length(swapworthy) == 1) {
      selected_precinct <- swapworthy[1]
    } else {
      selected_precinct <- sample(swapworthy, 1)
    }
    
    current_district <- current_data$district[selected_precinct]
    neighboring_districts_str <- current_data$neighboring_districts[selected_precinct]
    
    if (neighboring_districts_str == "") {
      next
    }
    
    target_districts <- as.numeric(strsplit(neighboring_districts_str, ",")[[1]])
    
    if (length(target_districts) == 0) {
      next
    }
    
    # Safe sampling for target
    if (length(target_districts) == 1) {
      target_district <- target_districts[1]
    } else {
      target_district <- sample(target_districts, 1)
    }
    
    # Execute swap
    current_data$district[selected_precinct] <- target_district
    swaps_made <- swaps_made + 1
    
    # Update analysis after each swap for next iteration
    if (i < actual_s) {
      current_data <- update_analysis(current_data, neighbours_queen, neighbours_rook)
    }
  }
  
  # Final analysis after all swaps
  current_data <- update_analysis(current_data, neighbours_queen, neighbours_rook)
  
  return(list(nc_data = current_data, swaps_made = swaps_made, actual_s = actual_s))
}

# Function: Create district map with dissolved boundaries
create_district_map <- function(nc_data, title) {
  # Create dissolved districts for boundaries
  districts_dissolved <- nc_data %>%
    group_by(district) %>%
    summarise(.groups = "drop")
  
  ggplot() +
    geom_sf(data = nc_data, aes(fill = factor(district)), color = "white", size = 0.02) +
    geom_sf(data = districts_dissolved, fill = NA, color = "black", size = 0.5) +
    scale_fill_discrete(name = "District") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = title)
}

# Global variable to store baseline metrics for percentage scaling
BASELINE_METRICS <- NULL

# Function: LIGHTNING-FAST district evaluation with PERCENTAGE-BASED SCALING + CONTIGUITY
evaluate_districts <- function(nc_data, 
                               neighbours_queen,
                               weights = c(1, 1, 10, 1, 5), # Added contiguity weight
                               lower_is_better = TRUE, 
                               verbose = FALSE,
                               return_metrics = FALSE) {
  
  # Validate inputs
  if (length(weights) != 5) {
    stop("weights must be a vector of length 5: [max_compactness, mean_compactness, max_pop_dev, mean_pop_dev, contiguity]")
  }
  
  required_cols <- c("total_p", "centroid_x", "centroid_y")
  missing_cols <- required_cols[!required_cols %in% names(nc_data)]
  if (length(missing_cols) > 0) {
    stop("nc_data must contain columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Get basic info
  districts <- unique(nc_data$district)
  n_districts <- length(districts)
  total_pop <- sum(nc_data$total_p, na.rm = TRUE)
  ideal_pop_per_district <- total_pop / n_districts
  
  # Pre-allocate vectors
  compactness_scores <- numeric(n_districts)
  district_pops <- numeric(n_districts)
  
  # Calculate metrics for each district using LIGHTNING-FAST bounding box method
  for (i in seq_along(districts)) {
    d <- districts[i]
    district_precincts <- nc_data$district == d
    
    if (sum(district_precincts) == 0) next
    
    # Population (simple sum)
    district_pops[i] <- sum(nc_data$total_p[district_precincts], na.rm = TRUE)
    
    # LIGHTNING-FAST bounding box compactness
    district_x_coords <- nc_data$centroid_x[district_precincts]
    district_y_coords <- nc_data$centroid_y[district_precincts]
    
    # Calculate bounding box dimensions
    length_x <- max(district_x_coords) - min(district_x_coords)
    length_y <- max(district_y_coords) - min(district_y_coords)
    
    # Handle edge case of single precinct districts
    if (length_x == 0) length_x <- 1
    if (length_y == 0) length_y <- 1
    
    # Aspect ratio: max dimension / min dimension
    # Perfect square = 1.0, elongated shapes > 1.0
    compactness_scores[i] <- max(length_x, length_y) / min(length_x, length_y)
  }
  
  # Calculate population deviations
  pop_deviations <- abs(district_pops - ideal_pop_per_district)
  
  # Calculate contiguity score
  contiguity_result <- calculate_contiguity_score(nc_data, neighbours_queen)
  contiguity_score <- contiguity_result$score
  
  # Compute the five key metrics
  max_compactness <- max(compactness_scores)      # Worst (least compact) district
  mean_compactness <- mean(compactness_scores)    # Average compactness
  max_pop_deviation <- max(pop_deviations)        # Worst population imbalance
  mean_pop_deviation <- mean(pop_deviations)      # Average population imbalance
  
  # Store raw metrics for progress tracking
  raw_metrics <- c(max_compactness, mean_compactness, max_pop_deviation, mean_pop_deviation, contiguity_score)
  names(raw_metrics) <- c("max_compactness", "mean_compactness", "max_pop_deviation", "mean_pop_deviation", "contiguity")
  
  # PERCENTAGE-BASED SCALING: Store baseline on first call
  if (is.null(BASELINE_METRICS)) {
    BASELINE_METRICS <<- raw_metrics
    if (verbose) {
      cat("🎯 BASELINE METRICS ESTABLISHED:\n")
      cat(sprintf("  Max Compactness: %.4f\n", BASELINE_METRICS["max_compactness"]))
      cat(sprintf("  Mean Compactness: %.4f\n", BASELINE_METRICS["mean_compactness"]))
      cat(sprintf("  Max Pop Deviation: %s\n", format(round(BASELINE_METRICS["max_pop_deviation"]), big.mark = ",")))
      cat(sprintf("  Mean Pop Deviation: %s\n", format(round(BASELINE_METRICS["mean_pop_deviation"]), big.mark = ",")))
      cat(sprintf("  Contiguity Score: %.4f (Components: %d/%d)\n", 
                  BASELINE_METRICS["contiguity"], contiguity_result$total_components, contiguity_result$target_components))
      cat("\n")
    }
  }
  
  # Calculate percentage-based penalties (1.0 = no improvement, 0.0 = perfect)
  # For compactness: What fraction of possible improvement remains?
  max_compactness_penalty <- pmax(0, (max_compactness - 1.0) / (BASELINE_METRICS["max_compactness"] - 1.0))
  mean_compactness_penalty <- pmax(0, (mean_compactness - 1.0) / (BASELINE_METRICS["mean_compactness"] - 1.0))
  
  # For population: What fraction of original deviation remains?
  max_pop_penalty <- max_pop_deviation / BASELINE_METRICS["max_pop_deviation"] 
  mean_pop_penalty <- mean_pop_deviation / BASELINE_METRICS["mean_pop_deviation"]
  
  # For contiguity: Convert to penalty (1.0 - score, since higher contiguity is better)
  contiguity_penalty <- 1.0 - contiguity_score
  
  # Apply power penalization (squaring for distance-like behavior)
  penalties <- c(
    max_compactness_penalty^2,
    mean_compactness_penalty^2,
    max_pop_penalty^2,
    mean_pop_penalty^2,
    contiguity_penalty^2
  )
  
  # Apply weights and sum
  weighted_penalties <- penalties * weights
  final_score <- sum(weighted_penalties)
  
  # Verbose output if requested
  if (verbose) {
    cat("=== PERCENTAGE-BASED DISTRICT EVALUATION + CONTIGUITY ===\n")
    cat("Districts:", n_districts, "\n")
    cat("Total population:", format(total_pop, big.mark = ","), "\n")
    cat("Ideal pop per district:", format(round(ideal_pop_per_district), big.mark = ","), "\n\n")
    
    cat("CONTIGUITY ANALYSIS:\n")
    cat(sprintf("  Connected components: %d (target: %d)\n", contiguity_result$total_components, contiguity_result$target_components))
    cat(sprintf("  Contiguity score: %.4f (%.1f%% of perfect contiguity)\n", contiguity_score, contiguity_score * 100))
    
    cat("\nBOUNDING BOX COMPACTNESS (1.0 = perfect square):\n")
    for (i in seq_along(districts)) {
      cat(sprintf("  District %d: %.4f (pop: %s)\n", 
                  districts[i], compactness_scores[i], 
                  format(district_pops[i], big.mark = ",")))
    }
    cat(sprintf("  Maximum (worst) compactness: %.4f\n", max_compactness))
    cat(sprintf("  Mean compactness: %.4f\n", mean_compactness))
    
    cat("\nPOPULATION DEVIATIONS:\n")
    for (i in seq_along(districts)) {
      deviation_pct <- (pop_deviations[i] / ideal_pop_per_district) * 100
      cat(sprintf("  District %d: %s (%.1f%% from ideal)\n", 
                  districts[i], format(round(pop_deviations[i]), big.mark = ","), deviation_pct))
    }
    cat(sprintf("  Maximum deviation: %s (%.1f%%)\n", 
                format(round(max_pop_deviation), big.mark = ","), 
                (max_pop_deviation / ideal_pop_per_district) * 100))
    cat(sprintf("  Mean deviation: %s (%.1f%%)\n", 
                format(round(mean_pop_deviation), big.mark = ","), 
                (mean_pop_deviation / ideal_pop_per_district) * 100))
    
    cat("\nPERCENTAGE-BASED PENALTY BREAKDOWN:\n")
    cat(sprintf("  Max compactness: %.1f%% remaining → penalty %.4f (weight=%.1f) → weighted=%.6f\n", 
                max_compactness_penalty*100, max_compactness_penalty^2, weights[1], weighted_penalties[1]))
    cat(sprintf("  Mean compactness: %.1f%% remaining → penalty %.4f (weight=%.1f) → weighted=%.6f\n", 
                mean_compactness_penalty*100, mean_compactness_penalty^2, weights[2], weighted_penalties[2]))
    cat(sprintf("  Max pop deviation: %.1f%% remaining → penalty %.4f (weight=%.1f) → weighted=%.6f\n", 
                max_pop_penalty*100, max_pop_penalty^2, weights[3], weighted_penalties[3]))
    cat(sprintf("  Mean pop deviation: %.1f%% remaining → penalty %.4f (weight=%.1f) → weighted=%.6f\n", 
                mean_pop_penalty*100, mean_pop_penalty^2, weights[4], weighted_penalties[4]))
    cat(sprintf("  Contiguity: %.1f%% penalty → penalty² %.4f (weight=%.1f) → weighted=%.6f\n", 
                contiguity_penalty*100, contiguity_penalty^2, weights[5], weighted_penalties[5]))
    cat(sprintf("  TOTAL SCORE: %.6f\n", final_score))
  }
  
  # Return score and optionally raw metrics
  if (return_metrics) {
    result <- list(
      score = if (lower_is_better) final_score else -final_score,
      metrics = raw_metrics,
      contiguity_components = contiguity_result$total_components,
      percentage_penalties = c(max_compactness_penalty, mean_compactness_penalty, max_pop_penalty, mean_pop_penalty, contiguity_penalty)
    )
    return(result)
  } else {
    return(if (lower_is_better) final_score else -final_score)
  }
}

# Helper function: Quick evaluation summary
evaluate_districts_summary <- function(nc_data, neighbours_queen, weights = c(1, 1, 10, 1, 5)) {
  result <- evaluate_districts(nc_data, neighbours_queen, weights = weights, verbose = TRUE, return_metrics = TRUE)
  return(result)
}

# Enhanced simulated annealing acceptance function
accept_worse_solution <- function(current_score, candidate_score, temperature, lower_is_better = TRUE) {
  if (lower_is_better) {
    delta <- candidate_score - current_score  # positive = worse
    if (delta <= 0) return(TRUE)  # Always accept improvements
    return(runif(1) < exp(-delta / temperature))  # Probabilistically accept worse
  } else {
    delta <- current_score - candidate_score  # positive = worse  
    if (delta <= 0) return(TRUE)
    return(runif(1) < exp(-delta / temperature))
  }
}

# Start  ---------------------------------------------------------------------

cat("Starting ENHANCED district optimization with CONTIGUITY BREAKING...\n\n")

# Reset baseline metrics for new run
BASELINE_METRICS <- NULL

nc <- initialize_districts_from_shapefile("geom/houstonia_precincts2.shp",
                                          d = 10
)

# Parameters - UPDATED with contiguity
s <- 5 # Maximum swaps per edit (replaces w)
r <- 2000 # Runs: number of edits
evaluation_weights <- c(5, 0.5, 25, 5, 15) # Weights for [max_compactness, mean_compactness, max_pop_dev, mean_pop_dev, contiguity]
lower_is_better <- TRUE # Whether lower scores are better

# Enhanced annealing parameters
initial_temp_factor <- 0.1  # Temperature as fraction of initial score
cooling_rate <- 0.9995      # Cool by 0.05% each iteration
min_temperature <- 1e-6     # Minimum temperature floor

# Get spatial neighbors
neighbours_queen <- sfdep::st_contiguity(nc, queen = TRUE)
neighbours_rook <- sfdep::st_contiguity(nc, queen = FALSE)

# Initialize base map and analysis
base_map <- nc
base_map <- update_analysis(base_map, neighbours_queen, neighbours_rook)

# Get detailed initial evaluation
cat("=== INITIAL DISTRICT CONFIGURATION ===\n")
initial_result <- evaluate_districts_summary(base_map, neighbours_queen, weights = evaluation_weights)
base_score <- initial_result$score
initial_metrics <- initial_result$metrics
initial_components <- initial_result$contiguity_components

cat("\nInitial comprehensive score:", round(base_score, 6), "\n")
cat("Border precincts:", sum(base_map$border), "\n")
cat("Connected components:", initial_components, "(target:", length(unique(base_map$district)), ")\n\n")

# Store initial map for comparison
initial_map <- create_district_map(base_map, "Initial Districts")

# Track progress for all 5 metrics
sample_interval <- 10  # Track every 10 runs
metrics_sample_runs <- seq(0, r, by = sample_interval)
metrics_history <- array(NA, dim = c(length(metrics_sample_runs), 5))
colnames(metrics_history) <- c("max_compactness", "mean_compactness", "max_pop_deviation", "mean_pop_deviation", "contiguity")

# Store initial metrics
initial_penalties <- initial_result$percentage_penalties
metrics_history[1, ] <- initial_penalties
metrics_counter <- 2

# Track progress
score_history <- base_score
improvements <- 0
total_swaps_made <- 0
temperature_history <- numeric(r)

# Calculate initial temperature
initial_temp <- base_score * initial_temp_factor

# Run Optimization --------------------------------------------------------

cat("=== STARTING ENHANCED OPTIMIZATION LOOP ===\n")
cat("Parameters:\n")
cat("  Max swaps per edit:", s, "\n")
cat("  Total runs:", r, "\n")
cat("  Evaluation weights: [max_compactness=", evaluation_weights[1], 
    ", mean_compactness=", evaluation_weights[2], 
    ", max_pop_dev=", evaluation_weights[3], " (HEAVILY WEIGHTED)", 
    ", mean_pop_dev=", evaluation_weights[4], 
    ", contiguity=", evaluation_weights[5], " (NEW)]\n")
cat("  Initial temperature:", round(initial_temp, 6), "\n")
cat("  Cooling rate:", cooling_rate, "\n")
cat("  Min temperature:", min_temperature, "\n\n")

cat("CONTIGUITY BREAKING ENABLED:\n")
cat("  ✅ Integral precincts filter REMOVED - all border precincts eligible\n")
cat("  ✅ Multiple swaps per edit (1 to", s, ") - breaking contiguity allowed\n")
cat("  ✅ Contiguity score integrated into optimization\n")
cat("  ✅ Enhanced simulated annealing with temperature tracking\n")
cat("  ✅ All 5 metrics tracked every", sample_interval, "runs for visualization\n\n")

for (run in 1:r) {
  # Calculate current temperature
  current_temp <- max(min_temperature, initial_temp * (cooling_rate ^ run))
  temperature_history[run] <- current_temp
  
  # Start from current base
  candidate_map <- base_map
  
  # Execute one edit (1 to s swaps randomly)
  edit_result <- execute_edit(candidate_map, neighbours_queen, neighbours_rook, s)
  candidate_map <- edit_result$nc_data
  swaps_in_run <- edit_result$swaps_made
  total_swaps_made <- total_swaps_made + swaps_in_run
  
  # Evaluate candidate configuration
  candidate_score <- evaluate_districts(candidate_map, neighbours_queen, weights = evaluation_weights, lower_is_better = lower_is_better)
  
  # Determine if improvement
  is_improvement <- if (lower_is_better) {
    candidate_score < base_score
  } else {
    candidate_score > base_score
  }
  
  # Accept improvements OR worse solutions probabilistically
  if (is_improvement || accept_worse_solution(base_score, candidate_score, current_temp, lower_is_better)) {
    base_map <- candidate_map
    base_score <- candidate_score
    if (is_improvement) improvements <- improvements + 1
  }
  
  score_history <- c(score_history, base_score)
  
  # Track metrics at sample intervals
  if (run %% sample_interval == 0 && metrics_counter <= nrow(metrics_history)) {
    current_result <- evaluate_districts(base_map, neighbours_queen, weights = evaluation_weights, 
                                         lower_is_better = lower_is_better, return_metrics = TRUE)
    metrics_history[metrics_counter, ] <- current_result$percentage_penalties
    metrics_counter <- metrics_counter + 1
  }
  
  # Enhanced progress reporting every 100 runs
  if (run %% 100 == 0) {
    cat("=== Run", run, "===\n")
    current_result <- evaluate_districts(base_map, neighbours_queen, weights = evaluation_weights, 
                                         lower_is_better = lower_is_better, return_metrics = TRUE)
    current_metrics <- current_result$metrics
    current_penalties <- current_result$percentage_penalties
    current_components <- current_result$contiguity_components
    
    # Calculate % improvements from baseline (negative = improvement)
    pct_improvements <- (1 - current_penalties[1:4]) * 100
    contiguity_pct <- current_metrics["contiguity"] * 100
    
    cat("TARGET MEASURES & PROGRESS:\n")
    cat(sprintf("  Max Compactness: %.4f (%.1f%% improved from baseline)\n", 
                current_metrics["max_compactness"], pct_improvements[1]))
    cat(sprintf("  Mean Compactness: %.4f (%.1f%% improved from baseline)\n", 
                current_metrics["mean_compactness"], pct_improvements[2]))
    cat(sprintf("  Max Pop Deviation: %s (%.1f%% improved) **HEAVILY WEIGHTED**\n", 
                format(round(current_metrics["max_pop_deviation"]), big.mark = ","), 
                pct_improvements[3]))
    cat(sprintf("  Mean Pop Deviation: %s (%.1f%% improved from baseline)\n", 
                format(round(current_metrics["mean_pop_deviation"]), big.mark = ","), 
                pct_improvements[4]))
    cat(sprintf("  Contiguity: %.4f (%.1f%% of perfect) Components: %d/%d **NEW METRIC**\n", 
                current_metrics["contiguity"], contiguity_pct, current_components, 
                length(unique(base_map$district))))
    cat(sprintf("  Overall Score: %.6f | Temp: %.2e | Improvements: %d/%d | Total Swaps: %d\n\n", 
                base_score, current_temp, improvements, run, total_swaps_made))
  }
}

# Final Evaluation and Visualization -------------------------------------

cat("=== FINAL DISTRICT CONFIGURATION ===\n")
final_result <- evaluate_districts_summary(base_map, neighbours_queen, weights = evaluation_weights)
final_metrics <- final_result$metrics
final_penalties <- final_result$percentage_penalties
final_components <- final_result$contiguity_components

# Calculate final % improvements from baseline
final_pct_improvements <- (1 - final_penalties[1:4]) * 100
final_contiguity_pct <- final_metrics["contiguity"] * 100

cat("\n=== OPTIMIZATION COMPLETE ===\n")
cat("FINAL TARGET MEASURES & IMPROVEMENTS:\n")
cat(sprintf("  Max Compactness: %.4f (%.1f%% improved from baseline)\n", 
            final_metrics["max_compactness"], final_pct_improvements[1]))
cat(sprintf("  Mean Compactness: %.4f (%.1f%% improved from baseline)\n", 
            final_metrics["mean_compactness"], final_pct_improvements[2]))
cat(sprintf("  Max Pop Deviation: %s (%.1f%% improved) **HEAVILY WEIGHTED**\n", 
            format(round(final_metrics["max_pop_deviation"]), big.mark = ","), 
            final_pct_improvements[3]))
cat(sprintf("  Mean Pop Deviation: %s (%.1f%% improved from baseline)\n", 
            format(round(final_metrics["mean_pop_deviation"]), big.mark = ","), 
            final_pct_improvements[4]))
cat(sprintf("  Contiguity: %.4f (%.1f%% of perfect) Components: %d/%d **CONTIGUITY BREAKING ENABLED**\n", 
            final_metrics["contiguity"], final_contiguity_pct, final_components, 
            length(unique(base_map$district))))

cat(sprintf("\nOverall Score: %.6f → %.6f\n", score_history[1], base_score))
cat("Total runs:", r, "\n")
cat("Total swaps attempted:", total_swaps_made, "\n")
cat("Improvements accepted:", improvements, "\n")
cat("Final temperature:", round(min(temperature_history), 8), "\n")

# Calculate overall improvement
if (lower_is_better) {
  improvement <- score_history[1] - base_score
  improvement_pct <- improvement / score_history[1] * 100
  cat(sprintf("Overall Improvement: %.6f (%.2f%% reduction)\n", improvement, improvement_pct))
} else {
  improvement <- base_score - score_history[1]
  improvement_pct <- improvement / score_history[1] * 100
  cat(sprintf("Overall Improvement: %.6f (%.2f%% increase)\n", improvement, improvement_pct))
}

# Create final map
final_map <- create_district_map(base_map, "Optimized Districts (Contiguity Breaking Enabled)")

# Display maps
print(initial_map)
print(final_map)

# Visualize Plots ----------------------------------------------------------

# Plot score progress
score_df <- data.frame(
  run = 0:r,
  score = score_history
)

progress_plot <- ggplot(score_df, aes(x = run, y = score)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 0.5) +
  theme_minimal() +
  labs(
    title = "Enhanced District Optimization Progress with Contiguity Breaking",
    subtitle = paste("Weights: [max_comp=", evaluation_weights[1], 
                     ", mean_comp=", evaluation_weights[2], 
                     ", max_pop_dev=", evaluation_weights[3], " (HEAVY)", 
                     ", mean_pop_dev=", evaluation_weights[4],
                     ", contiguity=", evaluation_weights[5], " (NEW)]"),
    x = "Run",
    y = "Comprehensive Score (Lower = Better)"
  )

print(progress_plot)

# Plot temperature schedule
temp_df <- data.frame(
  run = 1:r,
  temperature = temperature_history
)

temp_plot <- ggplot(temp_df, aes(x = run, y = temperature)) +
  geom_line(color = "red", size = 1) +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Simulated Annealing Temperature Schedule",
    subtitle = paste("Initial temp:", round(initial_temp, 6), "| Cooling rate:", cooling_rate, "| Min temp:", min_temperature),
    x = "Run",
    y = "Temperature (log scale)"
  )

print(temp_plot)

# All 5 Metrics Over Time Plot --------------------------------------------

cat("📊 Creating comprehensive metrics visualization...\n")

# Convert metrics history to plotting format
metrics_df <- data.frame(
  run = metrics_sample_runs[1:nrow(metrics_history)],
  metrics_history[1:nrow(metrics_history), ]
) %>%
  pivot_longer(cols = -run, names_to = "metric", values_to = "penalty") %>%
  mutate(
    metric = case_when(
      metric == "max_compactness" ~ "Max Compactness",
      metric == "mean_compactness" ~ "Mean Compactness", 
      metric == "max_pop_deviation" ~ "Max Population Deviation",
      metric == "mean_pop_deviation" ~ "Mean Population Deviation",
      metric == "contiguity" ~ "Contiguity (Inverted)",
      TRUE ~ metric
    ),
    improvement_score = 1 - penalty  # Convert to improvement (1 = perfect, 0 = baseline)
  )

# Create the comprehensive metrics plot
metrics_plot <- ggplot(metrics_df, aes(x = run, y = improvement_score, color = metric)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(
    "Max Compactness" = "#e41a1c",
    "Mean Compactness" = "#ff7f00", 
    "Max Population Deviation" = "#4daf4a",
    "Mean Population Deviation" = "#377eb8",
    "Contiguity (Inverted)" = "#984ea3"
  )) +
  scale_y_continuous(
    limits = c(-0.3, 1.1),
    breaks = seq(-0.2, 1.0, 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "District Optimization: All 5 Target Metrics Over Time",
    subtitle = "Higher = Better | 100% = Perfect | 0% = Baseline | Negative = Worse than Baseline",
    x = "Optimization Run",
    y = "Improvement Score (% of Perfect)",
    caption = paste("Tracked every", sample_interval, "runs |", 
                    "Pop Dev weight =", evaluation_weights[3], 
                    "| Contiguity weight =", evaluation_weights[5])
  )

print(metrics_plot)

# Create summary table of final improvements
final_improvements <- data.frame(
  Metric = c("Max Compactness", "Mean Compactness", "Max Population Deviation", 
             "Mean Population Deviation", "Contiguity"),
  Initial_Penalty = c(1.0, 1.0, 1.0, 1.0, 1.0),  # All start at baseline
  Final_Penalty = final_result$percentage_penalties,
  Improvement_Score = 1 - final_result$percentage_penalties,
  Improvement_Pct = (1 - final_result$percentage_penalties) * 100,
  Weight = evaluation_weights
)

cat("\n📊 FINAL METRICS SUMMARY TABLE:\n")
print(final_improvements, digits = 4, row.names = FALSE)

# Contiguity analysis summary
cat("\n🔗 CONTIGUITY ANALYSIS SUMMARY:\n")
initial_contiguity_result <- calculate_contiguity_score(nc, neighbours_queen)
cat(sprintf("  Initial: %.4f (%d components for %d districts)\n", 
            initial_contiguity_result$score, initial_contiguity_result$total_components, 
            length(unique(nc$district))))
cat(sprintf("  Final: %.4f (%d components for %d districts)\n", 
            final_metrics["contiguity"], final_components, 
            length(unique(base_map$district))))
contiguity_change <- final_metrics["contiguity"] - initial_contiguity_result$score
cat(sprintf("  Change: %.4f ", contiguity_change))
if (contiguity_change > 0) {
  cat("(IMPROVED contiguity)\n")
} else if (contiguity_change < 0) {
  cat("(WORSE contiguity - but may be acceptable trade-off)\n")
} else {
  cat("(NO CHANGE in contiguity)\n")
}
