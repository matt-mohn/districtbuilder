# Package Load ------------------------------------------------------------

library(sf)
library(sfdep)
library(ggplot2)
library(dplyr)
library(stringr)

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

# Function: Tarjan's algorithm for articulation points
find_articulation_points <- function(adj_list) {
  n <- length(adj_list)
  visited <- rep(FALSE, n)
  disc <- rep(0, n)
  low <- rep(0, n)
  parent <- rep(-1, n)
  ap <- rep(FALSE, n)
  time <- 0
  
  tarjan_dfs <- function(u) {
    children <- 0
    visited[u] <<- TRUE
    time <<- time + 1
    disc[u] <<- low[u] <<- time
    
    for (v in adj_list[[u]]) {
      if (length(v) == 0) next
      if (!visited[v]) {
        parent[v] <<- u
        children <- children + 1
        tarjan_dfs(v)
        low[u] <<- min(low[u], low[v])
        if (parent[u] == -1 && children > 1) {
          ap[u] <<- TRUE
        }
        if (parent[u] != -1 && low[v] >= disc[u]) {
          ap[u] <<- TRUE
        }
      } else if (v != parent[u]) {
        low[u] <<- min(low[u], disc[v])
      }
    }
  }
  
  for (i in 1:n) {
    if (!visited[i] && length(adj_list[[i]]) > 0) {
      tarjan_dfs(i)
    }
  }
  
  return(ap)
}

# Function: Calculate integral precincts (uses queen contiguity)
calculate_integral <- function(nc_data, neighbours_queen) {
  district_filtered_neighbours <- vector("list", nrow(nc_data))
  
  for (i in seq_along(neighbours_queen)) {
    current_district <- nc_data$district[i]
    same_district_neighbors <- neighbours_queen[[i]][nc_data$district[neighbours_queen[[i]]] == current_district]
    district_filtered_neighbours[[i]] <- same_district_neighbors
  }
  
  integral_logical <- find_articulation_points(district_filtered_neighbours)
  return(as.numeric(integral_logical))
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

# Function: Update integral and border for current district assignment
update_analysis <- function(nc_data, neighbours_queen, neighbours_rook) {
  nc_data$integral <- calculate_integral(nc_data, neighbours_queen)
  border_results <- calculate_border(nc_data, neighbours_rook)
  nc_data$border <- border_results$border
  nc_data$neighboring_districts <- border_results$neighboring_districts
  return(nc_data)
}

# Function: Execute a single swap
execute_single_swap <- function(nc_data) {
  swapworthy <- which(nc_data$border == 1 & nc_data$integral == 0)
  
  if (length(swapworthy) == 0) {
    return(list(nc_data = nc_data, swapped = FALSE))
  }
  
  # Safe sampling
  if (length(swapworthy) == 1) {
    selected_precinct <- swapworthy[1]
  } else {
    selected_precinct <- sample(swapworthy, 1)
  }
  
  current_district <- nc_data$district[selected_precinct]
  neighboring_districts_str <- nc_data$neighboring_districts[selected_precinct]
  
  if (neighboring_districts_str == "") {
    return(list(nc_data = nc_data, swapped = FALSE))
  }
  
  target_districts <- as.numeric(strsplit(neighboring_districts_str, ",")[[1]])
  
  if (length(target_districts) == 0) {
    return(list(nc_data = nc_data, swapped = FALSE))
  }
  
  # Safe sampling for target
  if (length(target_districts) == 1) {
    target_district <- target_districts[1]
  } else {
    target_district <- sample(target_districts, 1)
  }
  
  # Execute swap
  nc_data$district[selected_precinct] <- target_district
  
  return(list(nc_data = nc_data, swapped = TRUE))
}

# Function: Execute w swaps (one edit) - now with random width
execute_edit <- function(nc_data, neighbours_queen, neighbours_rook, max_w) {
  current_data <- nc_data
  
  # Randomly choose number of swaps between 1 and max_w (biased toward lower values)
  weights <- rev(2^(0:(max_w - 1))) # Creates weights favoring fewer swaps
  actual_w <- sample(1:max_w, 1, prob = weights)
  
  swaps_made <- 0
  
  for (i in 1:actual_w) {
    # Update analysis before each swap
    current_data <- update_analysis(current_data, neighbours_queen, neighbours_rook)
    
    # Try to swap
    swap_result <- execute_single_swap(current_data)
    current_data <- swap_result$nc_data
    
    if (swap_result$swapped) {
      swaps_made <- swaps_made + 1
    }
  }
  
  # Final analysis after all swaps
  current_data <- update_analysis(current_data, neighbours_queen, neighbours_rook)
  
  return(list(nc_data = current_data, swaps_made = swaps_made, actual_w = actual_w))
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

# Function: LIGHTNING-FAST district evaluation with PERCENTAGE-BASED SCALING
# Parameters:
#   nc_data: spatial data frame with district assignments, total_p, centroid_x, centroid_y columns
#   weights: vector of 4 weights for [max_compactness, mean_compactness, max_pop_dev, mean_pop_dev]
#   lower_is_better: logical, whether lower scores indicate better districts
#   verbose: logical, whether to print detailed breakdown
# Returns: numeric score for the configuration (lower = better by default)

evaluate_districts <- function(nc_data, 
                               weights = c(1, 1, 10, 1), 
                               lower_is_better = TRUE, 
                               verbose = FALSE,
                               return_metrics = FALSE) {
  
  # Validate inputs
  if (length(weights) != 4) {
    stop("weights must be a vector of length 4: [max_compactness, mean_compactness, max_pop_dev, mean_pop_dev]")
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
  
  # Compute the four key metrics
  max_compactness <- max(compactness_scores)      # Worst (least compact) district
  mean_compactness <- mean(compactness_scores)    # Average compactness
  max_pop_deviation <- max(pop_deviations)        # Worst population imbalance
  mean_pop_deviation <- mean(pop_deviations)      # Average population imbalance
  
  # Store raw metrics for progress tracking
  raw_metrics <- c(max_compactness, mean_compactness, max_pop_deviation, mean_pop_deviation)
  names(raw_metrics) <- c("max_compactness", "mean_compactness", "max_pop_deviation", "mean_pop_deviation")
  
  # PERCENTAGE-BASED SCALING: Store baseline on first call
  if (is.null(BASELINE_METRICS)) {
    BASELINE_METRICS <<- raw_metrics
    if (verbose) {
      cat("🎯 BASELINE METRICS ESTABLISHED:\n")
      cat(sprintf("  Max Compactness: %.4f\n", BASELINE_METRICS["max_compactness"]))
      cat(sprintf("  Mean Compactness: %.4f\n", BASELINE_METRICS["mean_compactness"]))
      cat(sprintf("  Max Pop Deviation: %s\n", format(round(BASELINE_METRICS["max_pop_deviation"]), big.mark = ",")))
      cat(sprintf("  Mean Pop Deviation: %s\n", format(round(BASELINE_METRICS["mean_pop_deviation"]), big.mark = ",")))
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
  
  # Apply power penalization (squaring for distance-like behavior)
  penalties <- c(
    max_compactness_penalty^2,
    mean_compactness_penalty^2,
    max_pop_penalty^2,
    mean_pop_penalty^2
  )
  
  # Apply weights and sum
  weighted_penalties <- penalties * weights
  final_score <- sum(weighted_penalties)
  
  # Verbose output if requested
  if (verbose) {
    cat("=== PERCENTAGE-BASED DISTRICT EVALUATION ===\n")
    cat("Districts:", n_districts, "\n")
    cat("Total population:", format(total_pop, big.mark = ","), "\n")
    cat("Ideal pop per district:", format(round(ideal_pop_per_district), big.mark = ","), "\n\n")
    
    cat("BOUNDING BOX COMPACTNESS (1.0 = perfect square):\n")
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
    cat(sprintf("  TOTAL SCORE: %.6f\n", final_score))
  }
  
  # Return score and optionally raw metrics
  if (return_metrics) {
    result <- list(
      score = if (lower_is_better) final_score else -final_score,
      metrics = raw_metrics,
      percentage_penalties = c(max_compactness_penalty, mean_compactness_penalty, max_pop_penalty, mean_pop_penalty)
    )
    return(result)
  } else {
    return(if (lower_is_better) final_score else -final_score)
  }
}

# Helper function: Quick evaluation summary
evaluate_districts_summary <- function(nc_data, weights = c(1, 1, 10, 1)) {
  result <- evaluate_districts(nc_data, weights = weights, verbose = TRUE, return_metrics = TRUE)
  return(result)
}

# Start  ---------------------------------------------------------------------

cat("Starting PERCENTAGE-BASED district optimization...\n\n")

# Reset baseline metrics for new run
BASELINE_METRICS <- NULL

nc <- initialize_districts_from_shapefile("geom/houstonia_precincts2.shp",
                                          d = 10
)

# Parameters - UPDATED
w <- 25 # Width: swaps per edit
r <- 1000 # Runs: number of edits
evaluation_weights <- c(5, 0.5, 25, 5) # Weights for [max_compactness, mean_compactness, max_pop_dev, mean_pop_dev] - HEAVILY WEIGHT POP DEVIATION
lower_is_better <- TRUE # Whether lower scores are better

# Get spatial neighbors
neighbours_queen <- sfdep::st_contiguity(nc, queen = TRUE)
neighbours_rook <- sfdep::st_contiguity(nc, queen = FALSE)

# Initialize base map and analysis
base_map <- nc
base_map <- update_analysis(base_map, neighbours_queen, neighbours_rook)

# Get detailed initial evaluation
cat("=== INITIAL DISTRICT CONFIGURATION ===\n")
initial_result <- evaluate_districts_summary(base_map, weights = evaluation_weights)
base_score <- initial_result$score
initial_metrics <- initial_result$metrics

cat("\nInitial comprehensive score:", round(base_score, 6), "\n")
cat("Integral precincts:", sum(base_map$integral), "\n")
cat("Border precincts:", sum(base_map$border), "\n")
cat("Swapworthy precincts:", sum(base_map$border == 1 & base_map$integral == 0), "\n\n")

# Store initial map for comparison
initial_map <- create_district_map(base_map, "Initial Districts")

# Track progress
score_history <- base_score
improvements <- 0
total_swaps_made <- 0

# Run Optimization --------------------------------------------------------

cat("=== STARTING OPTIMIZATION LOOP ===\n")
cat("Parameters:\n")
cat("  Max swaps per run:", w, "\n")
cat("  Total runs:", r, "\n")
cat("  Evaluation weights: [max_compactness=", evaluation_weights[1], 
    ", mean_compactness=", evaluation_weights[2], 
    ", max_pop_dev=", evaluation_weights[3], " (HEAVILY WEIGHTED)", 
    ", mean_pop_dev=", evaluation_weights[4], "]\n\n")

cat("PERCENTAGE-BASED SCALING CONFIRMED:\n")
cat("  ✅ All metrics start at 1.0 (100% remaining) and decrease toward 0.0 (perfect)\n")
cat("  ✅ Compactness: (current-1.0)/(baseline-1.0) = % of possible improvement remaining\n")
cat("  ✅ Population: current/baseline = % of original deviation remaining\n")
cat("  ✅ Perfect automatic scaling - both metrics on same 0.0-1.0 scale!\n\n")

for (run in 1:r) {
  # Start from current base
  candidate_map <- base_map
  
  # Execute one edit (1 to w swaps randomly)
  edit_result <- execute_edit(candidate_map, neighbours_queen, neighbours_rook, w)
  candidate_map <- edit_result$nc_data
  swaps_in_run <- edit_result$swaps_made
  total_swaps_made <- total_swaps_made + swaps_in_run
  
  # Evaluate candidate configuration
  candidate_score <- evaluate_districts(candidate_map, weights = evaluation_weights, lower_is_better = lower_is_better)
  
  # Add to your optimization loop
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
  
  # Replace your acceptance logic with:
  initial_temp <- base_score * 0.1  # Start temperature as 10% of initial score
  cooling_rate <- 0.995  # Cool by 0.5% each iteration
  temperature <- initial_temp * (cooling_rate ^ run)
  
  is_improvement <- if (lower_is_better) {
    candidate_score < base_score
  } else {
    candidate_score > base_score
  }
  
  # NEW: Accept improvements OR worse solutions probabilistically
  if (is_improvement || accept_worse_solution(base_score, candidate_score, temperature, lower_is_better)) {
    base_map <- candidate_map
    base_score <- candidate_score
    if (is_improvement) improvements <- improvements + 1
  }
  
  score_history <- c(score_history, base_score)
  
  # Progress reporting every 100 runs with percentage improvements
  if (run %% 100 == 0) {
    cat("=== Run", run, "===\n")
    current_result <- evaluate_districts(base_map, weights = evaluation_weights, 
                                         lower_is_better = lower_is_better, return_metrics = TRUE)
    current_metrics <- current_result$metrics
    current_penalties <- current_result$percentage_penalties
    
    # Calculate % improvements from baseline (negative = improvement)
    pct_improvements <- (1 - current_penalties) * 100
    
    cat("TARGET MEASURES & PROGRESS:\n")
    cat(sprintf("  Max Compactness: %.4f (%.1f%% improved from baseline) [worst district]\n", 
                current_metrics["max_compactness"], pct_improvements[1]))
    cat(sprintf("  Mean Compactness: %.4f (%.1f%% improved from baseline) [average district]\n", 
                current_metrics["mean_compactness"], pct_improvements[2]))
    cat(sprintf("  Max Pop Deviation: %s (%.1f%% improved from baseline) **HEAVILY WEIGHTED** [worst district]\n", 
                format(round(current_metrics["max_pop_deviation"]), big.mark = ","), 
                pct_improvements[3]))
    cat(sprintf("  Mean Pop Deviation: %s (%.1f%% improved from baseline) [average district]\n", 
                format(round(current_metrics["mean_pop_deviation"]), big.mark = ","), 
                pct_improvements[4]))
    cat(sprintf("  Overall Score: %.6f | Improvements: %d | Total Swaps: %d\n\n", 
                base_score, improvements, total_swaps_made))
  }
}

# Final Evaluation and Visualization -------------------------------------

cat("=== FINAL DISTRICT CONFIGURATION ===\n")
final_result <- evaluate_districts_summary(base_map, weights = evaluation_weights)
final_metrics <- final_result$metrics
final_penalties <- final_result$percentage_penalties

# Calculate final % improvements from baseline
final_pct_improvements <- (1 - final_penalties) * 100

cat("\n=== OPTIMIZATION COMPLETE ===\n")
cat("FINAL TARGET MEASURES & IMPROVEMENTS:\n")
cat(sprintf("  Max Compactness: %.4f (%.1f%% improved from baseline) [worst district]\n", 
            final_metrics["max_compactness"], final_pct_improvements[1]))
cat(sprintf("  Mean Compactness: %.4f (%.1f%% improved from baseline) [average district]\n", 
            final_metrics["mean_compactness"], final_pct_improvements[2]))
cat(sprintf("  Max Pop Deviation: %s (%.1f%% improved from baseline) **HEAVILY WEIGHTED** [worst district]\n", 
            format(round(final_metrics["max_pop_deviation"]), big.mark = ","), 
            final_pct_improvements[3]))
cat(sprintf("  Mean Pop Deviation: %s (%.1f%% improved from baseline) [average district]\n", 
            format(round(final_metrics["mean_pop_deviation"]), big.mark = ","), 
            final_pct_improvements[4]))

cat(sprintf("\nOverall Score: %.6f → %.6f\n", score_history[1], base_score))
cat("Total runs:", r, "\n")
cat("Total swaps attempted:", total_swaps_made, "\n")
cat("Improvements accepted:", improvements, "\n")

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
final_map <- create_district_map(base_map, "Optimized Districts")

# Display maps
print(initial_map)
print(final_map)

# Visualize Plot ----------------------------------------------------------

# Plot score progress
score_df <- data.frame(
  run = 0:r,
  score = score_history
)

progress_plot <- ggplot(score_df, aes(x = run, y = score)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  theme_minimal() +
  labs(
    title = "Percentage-Based District Optimization Progress",
    subtitle = paste("Weights: [max_comp=", evaluation_weights[1], 
                     ", mean_comp=", evaluation_weights[2], 
                     ", max_pop_dev=", evaluation_weights[3], " (HEAVY)", 
                     ", mean_pop_dev=", evaluation_weights[4], "] | All metrics 0.0-1.0 scale"),
    x = "Run",
    y = "Comprehensive Score (Lower = Better)"
  )

print(progress_plot)
