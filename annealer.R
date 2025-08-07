nc <- read_sf("geom/precincts.shp")  %>% mutate(district = as.integer(district)) %>%
  filter(district %in% c(24,33,30,32))

# Load required libraries
library(sf)
library(sfdep)
library(ggplot2)
library(dplyr)

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

# Function: Calculate total perimeter of dissolved districts
calculate_total_perimeter <- function(nc_data) {
  # Dissolve by district
  districts_dissolved <- nc_data %>%
    group_by(district) %>%
    summarise(.groups = 'drop')
  
  # Calculate perimeter using boundary length
  boundaries <- st_cast(districts_dissolved, "MULTILINESTRING")
  perimeters <- st_length(boundaries)
  total_perimeter <- sum(as.numeric(perimeters))
  
  return(total_perimeter)
}

# Function: Execute w swaps (one edit) - now with random width
execute_edit <- function(nc_data, neighbours_queen, neighbours_rook, max_w) {
  current_data <- nc_data
  
  # Randomly choose number of swaps between 1 and max_w (biased toward lower values)
  weights <- rev(2^(0:(max_w-1)))  # Creates weights: 1->8, 2->4, 3->2, 4->1 (for max_w=4)
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
    summarise(.groups = 'drop')
  
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

# MAIN OPTIMIZATION ALGORITHM
cat("Starting district optimization...\n\n")

# Parameters
w <- 5   # Width: swaps per edit
r <- 200  # Runs: number of edits

# Get spatial neighbors
neighbours_queen <- sfdep::st_contiguity(nc, queen = TRUE)
neighbours_rook <- sfdep::st_contiguity(nc, queen = FALSE)

# Initialize base map and analysis
base_map <- nc
base_map <- update_analysis(base_map, neighbours_queen, neighbours_rook)
base_perimeter <- calculate_total_perimeter(base_map)

cat("Initial state:\n")
cat("Total perimeter:", round(base_perimeter, 2), "\n")
cat("Integral precincts:", sum(base_map$integral), "\n")
cat("Border precincts:", sum(base_map$border), "\n")
cat("Swapworthy precincts:", sum(base_map$border == 1 & base_map$integral == 0), "\n\n")

# Store initial map for comparison
initial_map <- create_district_map(base_map, "Initial Districts")

# Track progress
perimeter_history <- base_perimeter
improvements <- 0
total_swaps_made <- 0

# Run optimization
for (run in 1:r) {
  cat("=== Run", run, "===\n")
  
  # Start from current base
  candidate_map <- base_map
  
  # Execute one edit (1 to w swaps randomly)
  edit_result <- execute_edit(candidate_map, neighbours_queen, neighbours_rook, w)
  candidate_map <- edit_result$nc_data
  swaps_in_run <- edit_result$swaps_made
  actual_w_used <- edit_result$actual_w
  total_swaps_made <- total_swaps_made + swaps_in_run
  
  # Calculate new perimeter
  candidate_perimeter <- calculate_total_perimeter(candidate_map)
  
  cat("Attempted swaps:", actual_w_used, "| Successful swaps:", swaps_in_run, "\n")
  cat("Candidate perimeter:", round(candidate_perimeter, 2), "\n")
  cat("Base perimeter:", round(base_perimeter, 2), "\n")
  
  # Check if improvement
  if (candidate_perimeter < base_perimeter) {
    # Accept improvement
    base_map <- candidate_map
    base_perimeter <- candidate_perimeter
    improvements <- improvements + 1
    cat("*** IMPROVEMENT ACCEPTED ***\n")
  } else {
    cat("No improvement - reverting to base\n")
  }
  
  perimeter_history <- c(perimeter_history, base_perimeter)
  cat("Current best perimeter:", round(base_perimeter, 2), "\n\n")
}

# Final results
final_map <- create_district_map(base_map, "Optimized Districts")

cat("=== OPTIMIZATION COMPLETE ===\n")
cat("Total runs:", r, "\n")
cat("Total swaps attempted:", total_swaps_made, "\n")
cat("Improvements accepted:", improvements, "\n")
cat("Initial perimeter:", round(perimeter_history[1], 2), "\n")
cat("Final perimeter:", round(base_perimeter, 2), "\n")
cat("Improvement:", round(perimeter_history[1] - base_perimeter, 2), 
    "(", round((perimeter_history[1] - base_perimeter) / perimeter_history[1] * 100, 2), "%)\n")

# Display maps
print(initial_map)
print(final_map)

# Plot perimeter progress
perimeter_df <- data.frame(
  run = 0:r,
  perimeter = perimeter_history
)

ggplot(perimeter_df, aes(x = run, y = perimeter)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  theme_minimal() +
  labs(
    title = "District Perimeter Optimization Progress",
    x = "Run",
    y = "Total Perimeter",
    subtitle = paste("w =", w, "swaps per run, r =", r, "runs")
  )