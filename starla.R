nc <- read_sf("geom/precincts.shp")  %>% mutate(district = as.integer(district)) %>%
  filter(district %in% c(7,9,38,9,18,29))

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

# Function: Calculate median isoperimetric quotient of districts
calculate_median_isoperimetric <- function(nc_data) {
  # Dissolve by district
  districts_dissolved <- nc_data %>%
    group_by(district) %>%
    summarise(.groups = 'drop')
  
  # Calculate isoperimetric quotient for each district
  areas <- as.numeric(st_area(districts_dissolved))
  boundaries <- st_cast(districts_dissolved, "MULTILINESTRING")
  perimeters <- as.numeric(st_length(boundaries))
  
  # Isoperimetric quotient = 4π * Area / Perimeter²
  # Perfect circle = 1, lower values = less compact
  iso_quotients <- (4 * pi * areas) / (perimeters^2)
  
  return(median(iso_quotients, na.rm = TRUE))
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

# MAIN OPTIMIZATION ALGORITHM - BEAM SEARCH
cat("Starting district optimization with beam search...\n\n")

# Parameters
w <- 10   # Width: max swaps per edit
c <- 7   # Children: number of children each parent generates
r <- 2000  # Total runs
beam_width <- 7  # Keep top 3 candidates each round

# Get spatial neighbors
neighbours_queen <- sfdep::st_contiguity(nc, queen = TRUE)
neighbours_rook <- sfdep::st_contiguity(nc, queen = FALSE)

# Initialize with base map
base_map <- nc
base_map <- update_analysis(base_map, neighbours_queen, neighbours_rook)
base_iso_quotient <- calculate_median_isoperimetric(base_map)

cat("Initial state:\n")
cat("Median isoperimetric quotient:", round(base_iso_quotient, 4), "\n")
cat("Integral precincts:", sum(base_map$integral), "\n")
cat("Border precincts:", sum(base_map$border), "\n")
cat("Swapworthy precincts:", sum(base_map$border == 1 & base_map$integral == 0), "\n\n")

# Store initial map
initial_map <- create_district_map(base_map, "Initial Districts")

# Initialize beam with base map
beam <- list(list(nc_data = base_map, iso_quotient = base_iso_quotient, generation = 0))
best_iso_history <- c(base_iso_quotient)
runs_completed <- 0

generation <- 1

while(runs_completed < r) {
  cat("=== Generation", generation, "===\n")
  cat("Current beam size:", length(beam), "\n")
  
  # Generate candidates from each beam member
  candidates <- list()
  
  for(i in seq_along(beam)) {
    parent <- beam[[i]]
    cat("Parent", i, "iso quotient:", round(parent$iso_quotient, 4), "\n")
    
    # Generate c children from each parent
    for(j in 1:c) {
      if(runs_completed >= r) break
      
      # Execute edit from this parent
      edit_result <- execute_edit(parent$nc_data, neighbours_queen, neighbours_rook, w)
      child_iso_quotient <- calculate_median_isoperimetric(edit_result$nc_data)
      
      candidates[[length(candidates) + 1]] <- list(
        nc_data = edit_result$nc_data,
        iso_quotient = child_iso_quotient,
        generation = generation,
        parent_id = i,
        swaps_made = edit_result$swaps_made,
        actual_w = edit_result$actual_w
      )
      
      runs_completed <- runs_completed + 1
      cat("  Child", j, "- Attempted:", edit_result$actual_w, "Successful:", edit_result$swaps_made, "Iso quotient:", round(child_iso_quotient, 4), "\n")
    }
    
    if(runs_completed >= r) break
  }
  
  # Add current beam members to candidates (so they can compete)
  for(beam_member in beam) {
    candidates[[length(candidates) + 1]] <- beam_member
  }
  
  # Sort all candidates by iso quotient (DESCENDING - higher is better)
  candidate_iso_quotients <- sapply(candidates, function(x) x$iso_quotient)
  top_indices <- order(candidate_iso_quotients, decreasing = TRUE)[1:min(beam_width, length(candidates))]
  beam <- candidates[top_indices]
  
  # Track best iso quotient
  best_current <- max(candidate_iso_quotients)
  best_iso_history <- c(best_iso_history, best_current)
  
  cat("Generation", generation, "complete. Best iso quotient:", round(best_current, 4), "\n")
  cat("Improvement from initial:", round(best_current - base_iso_quotient, 4), 
      "(", round((best_current - base_iso_quotient) / base_iso_quotient * 100, 2), "%)\n\n")
  
  generation <- generation + 1
}

# Get final best solution
best_solution <- beam[[1]]  # Already sorted, so first is best
final_map <- create_district_map(best_solution$nc_data, "Optimized Districts (Beam Search)")

cat("=== OPTIMIZATION COMPLETE ===\n")
cat("Total runs completed:", runs_completed, "\n")
cat("Generations:", generation - 1, "\n")
cat("Initial median iso quotient:", round(base_iso_quotient, 4), "\n")
cat("Final median iso quotient:", round(best_solution$iso_quotient, 4), "\n")
cat("Total improvement:", round(best_solution$iso_quotient - base_iso_quotient, 4), 
    "(", round((best_solution$iso_quotient - base_iso_quotient) / base_iso_quotient * 100, 2), "%)\n")

# Display maps
print(initial_map)
print(final_map)

# Plot progress
progress_df <- data.frame(
  generation = 0:(length(best_iso_history)-1),
  best_iso_quotient = best_iso_history
)

ggplot(progress_df, aes(x = generation, y = best_iso_quotient)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  theme_minimal() +
  labs(
    title = "District Compactness Optimization Progress (Beam Search)",
    x = "Generation",
    y = "Best Median Isoperimetric Quotient",
    subtitle = paste("Beam width =", beam_width, ", Children per parent =", c, ", Total runs =", runs_completed)
  )