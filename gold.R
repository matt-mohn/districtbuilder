nc <- read_sf("geom/precincts.shp")  %>% mutate(district = as.integer(district)) %>%
  filter(district %in% c(7, 38, 9, 18, 29))


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

# Function: Calculate integral precincts for current district assignment (uses queen contiguity)
calculate_integral <- function(nc_data, neighbours_queen) {
  # Filter neighbors to only include same-district connections
  district_filtered_neighbours <- vector("list", nrow(nc_data))
  
  for (i in seq_along(neighbours_queen)) {
    current_district <- nc_data$district[i]
    same_district_neighbors <- neighbours_queen[[i]][nc_data$district[neighbours_queen[[i]]] == current_district]
    district_filtered_neighbours[[i]] <- same_district_neighbors
  }
  
  # Find articulation points
  integral_logical <- find_articulation_points(district_filtered_neighbours)
  return(as.numeric(integral_logical))
}

# Function: Calculate border precincts and neighboring districts (uses rook contiguity)
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

# Function: Find and execute a single swap
execute_swap <- function(nc_data, neighbours_queen, neighbours_rook) {
  # Find swapworthy precincts (border==1, integral==0)
  swapworthy <- which(nc_data$border == 1 & nc_data$integral == 0)
  
  if (length(swapworthy) == 0) {
    cat("No swapworthy precincts found!\n")
    return(nc_data)
  }
  
  # Randomly select a swapworthy precinct (safe sampling)
  if (length(swapworthy) == 1) {
    selected_precinct <- swapworthy[1]
  } else {
    selected_precinct <- sample(swapworthy, 1)
  }
  
  # Get target districts from the neighboring_districts column
  current_district <- nc_data$district[selected_precinct]
  neighboring_districts_str <- nc_data$neighboring_districts[selected_precinct]
  
  if (neighboring_districts_str == "") {
    cat("Selected precinct has no neighboring districts!\n")
    return(nc_data)
  }
  
  # Parse the neighboring districts
  target_districts <- as.numeric(strsplit(neighboring_districts_str, ",")[[1]])
  
  if (length(target_districts) == 0) {
    cat("No target districts available for swap!\n")
    return(nc_data)
  }
  
  # Randomly select target district (safe sampling)
  if (length(target_districts) == 1) {
    target_district <- target_districts[1]
  } else {
    target_district <- sample(target_districts, 1)
  }
  
  # Execute the swap
  cat("Swapping precinct", selected_precinct, "from district", current_district, "to district", target_district, "\n")
  nc_data$district[selected_precinct] <- target_district
  
  return(nc_data)
}

# Function: Create district map
create_district_map <- function(nc_data, title) {
  ggplot(nc_data) +
    geom_sf(aes(fill = factor(district)), color = "white", size = 0.05) +
    scale_fill_discrete(name = "District") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = title)
}

# MAIN ANALYSIS
cat("Starting iterative swapping analysis...\n")

# Get spatial neighbors - two different contiguity definitions
neighbours_queen <- sfdep::st_contiguity(nc, queen = TRUE)   # For integral (includes corners)
neighbours_rook <- sfdep::st_contiguity(nc, queen = FALSE)   # For border (edges only)

# Store original district assignment
nc$original_district <- nc$district

# Initial analysis
nc <- update_analysis(nc, neighbours_queen, neighbours_rook)

# Create original map
original_map <- create_district_map(nc, "Original Districts")

cat("Initial state:\n")
cat("Integral precincts:", sum(nc$integral), "\n")
cat("Border precincts:", sum(nc$border), "\n")
cat("Swapworthy precincts (border=1, integral=0):", sum(nc$border == 1 & nc$integral == 0), "\n\n")

# Perform 20 swaps
for (swap_iteration in 1:100) {
  cat("=== Swap iteration", swap_iteration, "===\n")
  
  # Execute swap
  nc <- execute_swap(nc, neighbours_queen, neighbours_rook)
  
  # Update analysis
  nc <- update_analysis(nc, neighbours_queen, neighbours_rook)
  
  # Print status
  cat("After swap", swap_iteration, ":\n")
  cat("Integral precincts:", sum(nc$integral), "\n")
  cat("Border precincts:", sum(nc$border), "\n")
  cat("Swapworthy precincts:", sum(nc$border == 1 & nc$integral == 0), "\n\n")
}

# Create final map
final_map <- create_district_map(nc, "Districts After 20 Swaps")

# Display both maps
print(original_map)
print(final_map)

# Summary of changes
district_changes <- sum(nc$district != nc$original_district)
cat("Final summary:\n")
cat("Total precincts that changed districts:", district_changes, "\n")
cat("Percentage changed:", round(district_changes/nrow(nc) * 100, 2), "%\n")