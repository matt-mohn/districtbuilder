# Spatial Contiguity and Population Optimization using Simulated Annealing

# Load all packages -------------------------------------------------------

library(sf)
library(dplyr)
library(ggplot2)
library(viridis)
library(igraph)
library(gridExtra)
library(gganimate)
library(magick)
library(gifski)


# Configuration parameters ------------------------------------------------

### Topline settings
NUM_GROUPS <- 10
USE_POPULATION <- TRUE
set.seed(42)
shapefile.link <- "geom/houstonia_precincts3.shp"

### Annealing settings
min_temp <- 0.01
max_iterations <- 100000
start_temperature <- 10
end_temperature <- 0.05

start_cooling_rate <- (end_temperature / start_temperature)^(1/max_iterations)


### 0.01 at end, #1 at start,

### Dynamic annealing configurations
slow_cooling_rate = 0.9999
slow_cooling_trigger = 10000000 ## ignore for now

### Settings for border focus
max.border.focus.pct <- 0.75
min.border.focus.pct <- 0.01
midway.border.focus.pct <- 0.75

### Settings for island focus
max.island.focus.pct <- 0.25
min.island.focus.pct <- 0.01
midway.island.focus.pct <- 0.25

midway.hinge.pct <- 0.5

### Settings for allowing a 'dynamic' stop
kill.dynamic = FALSE
kill.dynamic.pct = 0.1

### Settings for allowing a 'static' stop
kill.static = FALSE
kill.static.num = 5000

### Settings for allowing a 'success' stop
kill.success = TRUE
kill.success.max_cohesive = 0.005
kill.success.max_population = 0.02
kill.success.max_external = 0.1

config.annealer.iofreq = 1000
config.chart.snapshotfreq = 10

### Art settings
MAKE_GIF <- TRUE
config.gif.fps = 2
config.gif.width = 14
config.gif.height = 12
config.gif.dpi = 100
config.gif.iofreq = 10
config.gif.margin = 10
config.gif.title.size = 14
config.gif.size = 0.1
config.gif.snapshotfreq = 2000
dissolved_gif = FALSE

map.counties = TRUE
continuity.surge = FALSE


# Helper functions (later for separate file) ------------------------------

# Fast function to calculate cohesive index using precomputed adjacency
calculate_cohesive_index_fast <- function(group_indices) {
  n_shapes_group <- length(group_indices)
  if (n_shapes_group <= 1) {
    return(0)
  } # Single shape or empty group is perfectly contiguous
  
  # Extract subgraph for this group
  subgraph <- induced_subgraph(adj_graph, group_indices)
  
  # Find connected components
  components <- components(subgraph)
  largest_component_size <- max(components$csize)
  
  # Calculate cohesive index
  cohesive_index <- (n_shapes_group - largest_component_size) / n_shapes_group
  return(cohesive_index)
}

# Function to calculate population index (PI) for a group
calculate_population_index <- function(group_indices) {
  if (length(group_indices) == 0) {
    return(1)
  } # Empty group gets maximum penalty
  
  group_population <- sum(shapes$population[group_indices], na.rm = TRUE)
  population_deviation <- abs(group_population - target_population_per_group) / target_population_per_group
  
  return(population_deviation)
}

# Fast function to calculate external interface index (XI) for a group
calculate_external_interface_index_fast <- function(group_indices, group_assignments) {
  if (length(group_indices) == 0) {
    return(1)
  } # Empty group gets maximum penalty
  
  # Vectorized approach: check which shapes have neighbors in different groups
  current_group <- group_assignments[group_indices[1]] # All shapes in group have same assignment
  
  # Count border shapes using vectorized operations
  border_count <- sum(sapply(group_indices, function(shape_idx) {
    neighbors <- neighbor_lists[[shape_idx]]
    if (length(neighbors) == 0) {
      return(FALSE)
    }
    any(group_assignments[neighbors] != current_group)
  }))
  
  # Calculate XI as proportion of shapes that border other groups
  xi <- border_count / length(group_indices)
  return(xi)
}

# Function to calculate overall score
calculate_score_fast <- function(group_assignments,continuity.surge = FALSE) {
  cohesive_indices <- numeric(NUM_GROUPS)
  population_indices <- numeric(NUM_GROUPS)
  external_interface_indices <- numeric(NUM_GROUPS)
  
  for (i in 1:NUM_GROUPS) {
    group_indices <- which(group_assignments == i)
    if (length(group_indices) > 0) {
      cohesive_indices[i] <- calculate_cohesive_index_fast(group_indices)
      population_indices[i] <- calculate_population_index(group_indices)
      external_interface_indices[i] <- calculate_external_interface_index_fast(group_indices, group_assignments)
    } else {
      cohesive_indices[i] <- 1 # Empty group gets penalized
      population_indices[i] <- 1 # Empty group gets penalized
      external_interface_indices[i] <- 1 # Empty group gets penalized
    }
  }
  
  avg_cohesive <- mean(cohesive_indices)
  min_cohesive <- min(cohesive_indices)
  max_cohesive <- max(cohesive_indices)
  
  avg_population <- mean(population_indices)
  min_population <- min(population_indices)
  max_population <- max(population_indices)
  
  avg_external <- mean(external_interface_indices)
  min_external <- min(external_interface_indices)
  max_external <- max(external_interface_indices)
  
    ci_score <- 0.1 * (avg_cohesive*100) + 0.9 * (max_cohesive*100)^2
    pi_score <- 0.1 * (avg_population*100) + 0.9 * (max_population*100)^2
    xi_score <- 0.1 * (avg_external*100) + 0.9 * (max_external*100)^2 
    
    
    if(continuity.surge==FALSE)
    {
      score <- 0.7 * ci_score + 0.2 * pi_score + 0.1 * xi_score
    }


    if(continuity.surge)
    {
      score <- 0.7 * ci_score + 0.2 * pi_score + 0.1 * xi_score
    }
    
    
  return(list(
    score = score,
    avg_cohesive = avg_cohesive,
    min_cohesive = min_cohesive,
    max_cohesive = max_cohesive,
    avg_population = avg_population,
    min_population = min_population,
    max_population = max_population,
    avg_external = avg_external,
    min_external = min_external,
    max_external = max_external,
    cohesive_indices = cohesive_indices,
    population_indices = population_indices,
    external_interface_indices = external_interface_indices
  ))
}





# Run ---------------------------------------------------------------------

cat("Reading shapefile...\n")
# Read the shapefile
shapes <- st_read(shapefile.link, quiet = TRUE)


if("countyid"%in%names(shapes) && map.counties){
  
  counties.shape <- shapes %>% 
    group_by(countyid) %>% 
    summarise(geometry = st_union(geometry))
  counties <- geom_sf(data=counties.shape,
                      fill=NA,
                      color="white",
                      linewidth=0.5)
  
}else{
  map.counties <- FALSE
  counties <- NULL
}

# Parse population data
cat("Parsing population data...\n")
shapes$population <- as.integer(shapes$total_p)
total_population <- sum(shapes$population, na.rm = TRUE)
target_population_per_group <- total_population / NUM_GROUPS

cat("Total population:", total_population, "\n")
cat("Target population per group:", round(target_population_per_group), "\n")
cat("Number of groups:", NUM_GROUPS, "\n")
cat("Using population optimization:", USE_POPULATION, "\n\n")

cat("Computing spatial adjacency matrix (this may take a moment)...\n")
# Compute adjacency matrix ONCE at the beginning
# Calculate shared boundary length and only count as adjacent if > some threshold
# Only edges that share a line segment (not just corner touches)
adj_matrix <- st_relate(shapes, shapes, pattern = "F***1****", sparse = FALSE)
diag(adj_matrix) <- FALSE

n_shapes <- nrow(shapes)

cat("Creating adjacency graph...\n")
# Create igraph object for fast connected components calculation
adj_graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# Precompute neighbor lists for fast XI calculation
neighbor_lists <- vector("list", n_shapes)
for (i in 1:n_shapes) {
  neighbor_lists[[i]] <- which(adj_matrix[i, ])
}


cat("Initializing random group assignments...\n")
# Initialize random group assignments
current_groups <- sample(1:NUM_GROUPS, n_shapes, replace = TRUE)

# Ensure each group has at least one shape
for (i in 1:NUM_GROUPS) {
  if (sum(current_groups == i) == 0) {
    random_shape <- sample(1:n_shapes, 1)
    current_groups[random_shape] <- i
  }
}

# Calculate initial score
current_metrics <- calculate_score_fast(current_groups)
current_score <- current_metrics$score

cat("Initial Score:", round(current_score, 4), "\n")
cat("Initial Avg Cohesive Index:", round(current_metrics$avg_cohesive, 4), "\n")
cat("Initial Min Cohesive Index:", round(current_metrics$min_cohesive, 4), "\n")
cat("Initial Max Cohesive Index:", round(current_metrics$max_cohesive, 4), "\n")
cat("Initial Avg External Interface Index:", round(current_metrics$avg_external, 4), "\n")
cat("Initial Min External Interface Index:", round(current_metrics$min_external, 4), "\n")
cat("Initial Max External Interface Index:", round(current_metrics$max_external, 4), "\n")
  cat("Initial Avg Population Index:", round(current_metrics$avg_population, 4), "\n")
  cat("Initial Min Population Index:", round(current_metrics$min_population, 4), "\n")
  cat("Initial Max Population Index:", round(current_metrics$max_population, 4), "\n")
cat("\n")

# Store initial state for plotting
initial_groups <- current_groups

initial_temp <- start_temperature
cooling_rate <- start_cooling_rate

cat("Using dynamic parameters:\n")
cat("Max iterations:", max_iterations, "\n")
cat("Initial temperature:", initial_temp, "\n")
cat("Cooling rate:", cooling_rate, "\n\n")

# Storage for tracking progress
progress_data <- data.frame(
  iteration = integer(),
  avg_cohesive = numeric(),
  min_cohesive = numeric(),
  max_cohesive = numeric(),
  avg_population = numeric(),
  min_population = numeric(),
  max_population = numeric(),
  avg_external = numeric(),
  min_external = numeric(),
  max_external = numeric(),
  score = numeric(),
  temperature = numeric()
)

# Storage for gif snapshots if enabled
if (MAKE_GIF) {
  gif_snapshots <- data.frame(
    iteration = integer(),
    shape_id = integer(),
    group = integer()
  )
  snapshot_counter <- 0
}

# Simulated annealing loop
temperature <- initial_temp
best_groups <- current_groups
best_score <- current_score
best_metrics <- current_metrics
no_improvement_count <- 0

cat("Starting simulated annealing...\n")
start_time <- Sys.time()

continuity.surge <- F

for (iter in 1:max_iterations) {
  # Calculate probability for border-focused and island-focused edits
  midway_point <- as.integer(round(max_iterations*midway.hinge.pct,digits=0))
  p_border <- min(max.border.focus.pct, min.border.focus.pct + ((midway.border.focus.pct - min.border.focus.pct) * min(iter / midway_point, 1)))
  p_island <- min(max.island.focus.pct , min.island.focus.pct + ((midway.island.focus.pct - min.island.focus.pct) * min(iter / midway_point, 1)))
  
  if(continuity.surge)
  {
    p_border <- p_border^(1/3)
    p_island <- p_island^(1/3)
  }
  
  
  # Determine edit type
  rand_val <- runif(1)
  use_border_focus <- rand_val < p_border
  use_island_focus <- use_border_focus && runif(1) < p_island
  
  # Single random swap per iteration with potential focus
  candidate_groups <- current_groups
  
  if (use_island_focus) {
    # Island-focused edit: find shapes that don't border any shapes from their own group
    island_shapes <- c()
    for (shape_idx in 1:n_shapes) {
      current_group <- candidate_groups[shape_idx]
      neighbors <- neighbor_lists[[shape_idx]]
      if (length(neighbors) > 0) {
        same_group_neighbors <- sum(candidate_groups[neighbors] == current_group)
        if (same_group_neighbors == 0) { # No neighbors from same group = island
          island_shapes <- c(island_shapes, shape_idx)
        }
      } else {
        # Shape with no neighbors is also considered an island
        island_shapes <- c(island_shapes, shape_idx)
      }
    }
    
    if (length(island_shapes) > 0) {
      # Select random island shape
      shape_idx <- if (length(island_shapes) == 1) island_shapes[1] else sample(island_shapes, 1)
      current_group <- candidate_groups[shape_idx]
      
      # Find which groups this shape borders
      neighbors <- neighbor_lists[[shape_idx]]
      if (length(neighbors) > 0) {
        neighbor_groups <- unique(candidate_groups[neighbors])
        available_groups <- setdiff(neighbor_groups, current_group)
      } else {
        available_groups <- setdiff(1:NUM_GROUPS, current_group)
      }
      
      if (length(available_groups) > 0) {
        new_group <- if (length(available_groups) == 1) available_groups[1] else sample(available_groups, 1)
        candidate_groups[shape_idx] <- new_group
      } else {
        # Fall back to regular random edit if no available groups
        use_border_focus <- FALSE
        use_island_focus <- FALSE
      }
    } else {
      # No islands found, fall back to border-focused
      use_island_focus <- FALSE
    }
  }
  
  if (use_border_focus && !use_island_focus) {
    # Border-focused edit: find shapes that border other groups
    border_shapes <- c()
    for (shape_idx in 1:n_shapes) {
      current_group <- candidate_groups[shape_idx]
      neighbors <- neighbor_lists[[shape_idx]]
      if (length(neighbors) > 0) {
        neighbor_groups <- candidate_groups[neighbors]
        if (any(neighbor_groups != current_group)) {
          border_shapes <- c(border_shapes, shape_idx)
        }
      }
    }
    
    if (length(border_shapes) > 0) {
      # Select random border shape
      shape_idx <- if (length(border_shapes) == 1) border_shapes[1] else sample(border_shapes, 1)
      current_group <- candidate_groups[shape_idx]
      
      # Find which groups this shape borders (excluding its current group)
      neighbors <- neighbor_lists[[shape_idx]]
      neighbor_groups <- unique(candidate_groups[neighbors])
      available_groups <- setdiff(neighbor_groups, current_group)
      
      if (length(available_groups) > 0) {
        new_group <- if (length(available_groups) == 1) available_groups[1] else sample(available_groups, 1)
        candidate_groups[shape_idx] <- new_group
      } else {
        # Fall back to regular random edit if no available groups
        use_border_focus <- FALSE
      }
    } else {
      # No border shapes found, fall back to regular random edit
      use_border_focus <- FALSE
    }
  }
  
  if (!use_border_focus) {
    # Regular random edit
    shape_idx <- sample(1:n_shapes, 1)
    current_group <- candidate_groups[shape_idx]
    available_groups <- setdiff(1:NUM_GROUPS, current_group)
    new_group <- if (length(available_groups) == 1) available_groups[1] else sample(available_groups, 1)
    candidate_groups[shape_idx] <- new_group
  }
  
  last.reading <- continuity.surge
  
  candidate_metrics <- calculate_score_fast(candidate_groups,continuity.surge)
  candidate_score <- candidate_metrics$score

  # Acceptance criteria
  delta <- candidate_score - current_score
  
  pass.through <- T
  
  if(last.reading) ## if contiguity mode activated
  {
    if(candidate_metrics$max_cohesive > current_metrics$max_cohesive)
    {
      pass.through <- F
    }
  }
  
  if (pass.through && (delta < 0 || (temperature > min_temp && runif(1) < exp(-delta / temperature)))) {
    # Accept the candidate solution
    current_groups <- candidate_groups
    current_score <- candidate_score
    current_metrics <- candidate_metrics
    
    # Update best solution if this is better
    if (current_score < best_score) {
      best_groups <- current_groups
      best_score <- current_score
      best_metrics <- current_metrics
      no_improvement_count <- 0
      
      ### If currently set OFF, and falls into range
      if(!last.reading && current_metrics$max_cohesive<0.3 && current_metrics$max_cohesive>0.03)
      {
        continuity.surge = TRUE
        cat("Entered contiguity mode")
        print(iter)
      }
      
      ### If currently set ON, and falls out of range
      if(last.reading && current_metrics$max_cohesive<0.03)
      {
        continuity.surge = FALSE
        cat("Left contiguity mode")
        print(iter)
      }
      
    } else {
      no_improvement_count <- no_improvement_count + 1
    }
  } else {
    no_improvement_count <- no_improvement_count + 1
  }
  
  # Adaptive cooling - slow down cooling if no improvement
  if (no_improvement_count > slow_cooling_trigger) {
    temperature <- temperature * slow_cooling_rate
  } else {
    temperature <- temperature * cooling_rate # Normal cooling
  }
  
  # Store progress every few iterations
  if (iter %% config.chart.snapshotfreq == 0) {
    progress_data <- rbind(progress_data, data.frame(
      iteration = iter,
      avg_cohesive = current_metrics$avg_cohesive,
      min_cohesive = current_metrics$min_cohesive,
      max_cohesive = current_metrics$max_cohesive,
      avg_population = current_metrics$avg_population,
      min_population = current_metrics$min_population,
      max_population = current_metrics$max_population,
      avg_external = current_metrics$avg_external,
      min_external = current_metrics$min_external,
      max_external = current_metrics$max_external,
      score = current_score,
      temperature = temperature
    ))
    
    # Take gif snapshot every few iterations
    if (MAKE_GIF && iter %% config.gif.snapshotfreq == 0) {
      snapshot_counter <- snapshot_counter + 1
      snapshot_data <- data.frame(
        iteration = iter,
        shape_id = 1:n_shapes,
        group = current_groups
      )
      gif_snapshots <- rbind(gif_snapshots, snapshot_data)
      cat(sprintf("Snapshot %d taken at iteration %d\n", snapshot_counter, iter))
    }
    
    if (iter %% config.annealer.iofreq == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      improvement <- round((progress_data$score[1] - current_score) / progress_data$score[1] * 100, 1)
      

        cat(sprintf(
          "Iter %d: Score=%.4f (%.1f%% better), CI: Avg=%.3f Max=%.3f, PI: Avg=%.3f Max=%.3f, XI: Avg=%.3f Max=%.3f, T=%.4f (%.1fs)\n",
          iter, current_score, improvement,
          current_metrics$avg_cohesive, current_metrics$max_cohesive,
          current_metrics$avg_population, current_metrics$max_population,
          current_metrics$avg_external, current_metrics$max_external,
          temperature, elapsed
        ))
       
    }
  }
  
  # Check for excellent solutions
  good_contiguity <- current_metrics$max_cohesive < kill.success.max_cohesive
  good_population <- current_metrics$max_population < kill.success.max_population
  good_compactness <- current_metrics$max_external < kill.success.max_external
  
  if (kill.success && good_contiguity && good_population && good_compactness) {
    cat(sprintf("Successful solution reached at iteration %d!\n", iter))
    cat(sprintf(
      "Max CI = %.4f, Max PI = %.4f, Max XI = %.4f\n",
      current_metrics$max_cohesive, current_metrics$max_population, current_metrics$max_external
    ))
    
    ### Add a gif if we haven't already
    if (MAKE_GIF && iter %% config.gif.snapshotfreq != 0) {
      snapshot_counter <- snapshot_counter + 1
      snapshot_data <- data.frame(
        iteration = iter,
        shape_id = 1:n_shapes,
        group = current_groups
      )
      gif_snapshots <- rbind(gif_snapshots, snapshot_data)
      cat(sprintf("Snapshot %d taken at iteration %d\n", snapshot_counter, iter))
    }
    
    
    
    
    break
  }
  
  # Early stopping if no improvement for a long time
  if (kill.static && no_improvement_count > kill.static.num) {
    cat(sprintf("No improvement for %d iterations, static count exceeded, stopping early at iteration %d\n", no_improvement_count, iter))
    break
  }
  
  # Early stopping if no improvement for a long time
  if (kill.dynamic && no_improvement_count > kill.dynamic.pct*max_iterations) {
    cat(sprintf("No improvement for %d iterations, dynamic count exceeded, stopping early at iteration %d\n", no_improvement_count, iter))
    break
  }
  
  
  
}

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("\nOptimization completed in", round(total_time, 1), "seconds\n")
cat("\nFinal Results:\n")
cat("Best Score:", round(best_score, 4), "\n")
cat("Best Avg Cohesive Index:", round(best_metrics$avg_cohesive, 4), "\n")
cat("Best Min Cohesive Index:", round(best_metrics$min_cohesive, 4), "\n")
cat("Best Max Cohesive Index:", round(best_metrics$max_cohesive, 4), "\n")
cat("Best Avg External Interface Index:", round(best_metrics$avg_external, 4), "\n")
cat("Best Min External Interface Index:", round(best_metrics$min_external, 4), "\n")
cat("Best Max External Interface Index:", round(best_metrics$max_external, 4), "\n")
  cat("Best Avg Population Index:", round(best_metrics$avg_population, 4), "\n")
  cat("Best Min Population Index:", round(best_metrics$min_population, 4), "\n")
  cat("Best Max Population Index:", round(best_metrics$max_population, 4), "\n")


improvement_pct <- round((progress_data$score[1] - best_score) / progress_data$score[1] * 100, 1)
cat("Overall improvement:", improvement_pct, "%\n")

# Show group sizes, populations, and indices
cat("\nFinal group statistics:\n")
for (i in 1:NUM_GROUPS) {
  group_indices <- which(best_groups == i)
  group_size <- length(group_indices)
  group_population <- sum(shapes$population[group_indices], na.rm = TRUE)
  group_cohesive <- best_metrics$cohesive_indices[i]
  group_pop_index <- best_metrics$population_indices[i]
  group_external <- best_metrics$external_interface_indices[i]
  pop_pct_target <- round(group_population / target_population_per_group * 100, 1)
  
  cat(sprintf(
    "Group %d: %d shapes, %d people (%.1f%% of target), CI=%.3f, PI=%.3f, XI=%.3f\n",
    i, group_size, group_population, pop_pct_target, group_cohesive, group_pop_index, group_external
  ))
}

# Prepare data for plotting
shapes$initial_group <- as.factor(initial_groups)
shapes$final_group <- as.factor(best_groups)

cat("\nGenerating plots...\n")

# Plot 1: Initial map
p1 <- ggplot(shapes) +
  geom_sf(aes(fill = initial_group), color = NA, size = 0.1) +
  scale_fill_viridis_d(name = "Group") +
  counties +
  ggtitle("Initial Random Assignment") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  )

# Plot 2: Final optimized map
p2 <- ggplot(shapes) +
  geom_sf(aes(fill = final_group), color = NA, size = 0.1) +
  scale_fill_viridis_d(name = "Group") +
  counties+
  ggtitle(paste0("Final Optimized Assignment (", improvement_pct, "% better)")) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  )

# Plot 3: Cohesive Index Progress
p3 <- ggplot(progress_data) +
  geom_line(aes(x = iteration, y = avg_cohesive, color = "CI Average"), size = 0.8) +
  geom_line(aes(x = iteration, y = min_cohesive, color = "CI Minimum"), size = 0.8) +
  geom_line(aes(x = iteration, y = max_cohesive, color = "CI Maximum"), size = 0.8) +
  scale_color_manual(values = c("CI Average" = "blue", "CI Minimum" = "green", "CI Maximum" = "red")) +
  labs(
    x = "Iteration",
    y = "Cohesive Index",
    title = "Cohesive Index Progress",
    color = "Metric"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  )

# Plot 4: Population Index Progress
p4 <- ggplot(progress_data) +
  geom_line(aes(x = iteration, y = avg_population, color = "PI Average"), size = 0.8) +
  geom_line(aes(x = iteration, y = min_population, color = "PI Minimum"), size = 0.8) +
  geom_line(aes(x = iteration, y = max_population, color = "PI Maximum"), size = 0.8) +
  scale_color_manual(values = c("PI Average" = "purple", "PI Minimum" = "orange", "PI Maximum" = "darkred")) +
  labs(
    x = "Iteration",
    y = "Population Index",
    title = "Population Index Progress",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  )

# Plot 5: External Interface Index Progress
p5 <- ggplot(progress_data) +
  geom_line(aes(x = iteration, y = avg_external, color = "XI Average"), size = 0.8) +
  geom_line(aes(x = iteration, y = min_external, color = "XI Minimum"), size = 0.8) +
  geom_line(aes(x = iteration, y = max_external, color = "XI Maximum"), size = 0.8) +
  scale_color_manual(values = c("XI Average" = "brown", "XI Minimum" = "pink", "XI Maximum" = "darkgreen")) +
  labs(
    x = "Iteration",
    y = "External Interface Index",
    title = "External Interface Index Progress",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  )

# Plot 6: Temperature progress
p6 <- ggplot(progress_data) +
  geom_line(aes(x = iteration, y = temperature), color = "gray50", size = 0.8) +
  labs(
    x = "Iteration",
    y = "Temperature",
    title = "Temperature Progress"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12)
  )


  combined_plot <- grid.arrange(p1, p2, p3+ylim(0,1)+
                                  geom_hline(yintercept=0)+
                                  geom_hline(yintercept=1), p4+ylim(0,1)+
                                  geom_hline(yintercept=0)+
                                  geom_hline(yintercept=1), p5+ylim(0,1)+
                                  geom_hline(yintercept=0)+
                                  geom_hline(yintercept=1), p6+ylim(0,start_temperature)+
                                  geom_hline(yintercept=start_temperature)+
                                  geom_hline(yintercept=0), ncol = 2, nrow = 3)


print(combined_plot)

# Create animated GIF if requested
if (MAKE_GIF && nrow(gif_snapshots) > 0) {
  cat("\nCreating animated GIF...\n")
  
  # Add final snapshot
  final_snapshot <- data.frame(
    iteration = max(progress_data$iteration),
    shape_id = 1:n_shapes,
    group = best_groups
  )
  gif_snapshots <- rbind(gif_snapshots, final_snapshot)
  
  # Get unique iterations for creating individual plots
  unique_iterations <- sort(unique(gif_snapshots$iteration))
  cat("Creating", length(unique_iterations), "frames...\n")
  
  # Create individual plots for each iteration (much faster than gganimate)
  plot_list <- list()
  
  for (i in seq_along(unique_iterations)) {
    iter <- unique_iterations[i]
    
    # Get group assignments for this iteration
    iter_data <- gif_snapshots[gif_snapshots$iteration == iter, ]
    
    # Create a temporary shapes object with current groups
    temp_shapes <- shapes
    temp_shapes$current_group <- as.factor(iter_data$group[match(1:n_shapes, iter_data$shape_id)])
    
    # Create plot for this iteration
    if(dissolved_gif == FALSE){
      p <- ggplot(temp_shapes) +
        geom_sf(aes(fill = as.factor(current_group)), color = NA, size = config.gif.size) +
        scale_fill_brewer(name = "district",palette = "Set3") +
        counties+
        ggtitle(paste("Optimization Progress: Iteration", iter)) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = config.gif.title.size, face = "bold"),
          legend.position = "bottom",
          plot.margin = margin(config.gif.margin, config.gif.margin, config.gif.margin, config.gif.margin)
        )
    }
    
    # Save plot as temporary PNG
    temp_file <- paste0("temp_frame_", sprintf("%03d", i), ".png")
    ggsave(temp_file, p, width = config.gif.width , 
           height = config.gif.height,
           dpi = config.gif.dpi)
    plot_list[[i]] <- temp_file
    
    if (i %% config.gif.iofreq == 0) cat("Created frame", i, "of", length(unique_iterations), "\n")
  }
  
  # Use magick to create GIF from individual images
  cat("Assembling GIF from frames...\n")
  imgs <- image_read(unlist(plot_list))
  gif <- image_animate(imgs, fps = config.gif.fps)
  image_write(gif, "optimization_progress.gif")
  
  # Clean up temporary files
  file.remove(unlist(plot_list))
  
  cat("Animated GIF saved as 'optimization_progress.gif'\n")
  cat("Total frames:", length(unique_iterations), "\n")
  cat("File size:", round(file.size("optimization_progress.gif") / 1024 / 1024, 1), "MB\n")
}

# Summary statistics
cat("\nSummary Statistics:\n")
cat("Total iterations run:", max(progress_data$iteration), "\n")
cat("Runtime:", round(total_time, 1), "seconds\n")
cat("Avg iterations per second:", round(max(progress_data$iteration) / total_time, 0), "\n")
cat("Mode: ","Contiguity + Population + Compactness", "\n")
