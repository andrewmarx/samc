# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

# This script is the source of the code for maze example vignettes
# The lines with `@knitr` are for processing by knitr

# TODO Incorporate more comments to make the script more useful on it's own
# without the vignette


## @knitr Part1
#
# Part 1 ----
#

## @knitr 1_library_1
library(samc)
library(raster)
library(gdistance)
library(viridisLite)


## @knitr 1_setup_1
maze_plot <- function(map, title, colors) {
  # start = 1 (top left), finish = last element (bottom right)
  sf <- xyFromCell(map, c(1, length(map)))

  plot(map, main = title, col = colors, axes = FALSE, box = FALSE, asp = 1)
  plot(rasterToPolygons(map), border = 'black', lwd = 1, add = TRUE)
  points(sf, pch = c('S', 'F'), cex = 1, font = 2)
}


## @knitr 1_setup_2
# A simple color palette with 2 colors
vir_col <- viridis(3)[2:3]


## @knitr 1_setup_3
maze = matrix(
  c(1,0,0,0,0,1,1,1,1,1,1,0,1,0,0,0,0,1,0,1,
    1,1,1,1,0,1,0,1,0,0,1,1,1,0,1,1,0,1,1,1,
    1,0,0,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0,
    1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,
    0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,1,
    1,1,0,1,0,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,
    0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,0,1,0,1,
    1,1,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,1,1,1,
    0,1,0,1,0,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,
    1,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,1,
    0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,
    1,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,1,
    0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,
    0,0,0,0,1,0,0,1,1,1,0,1,0,1,1,1,0,0,1,1,
    1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,1,0,
    1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,1,
    1,1,1,1,0,1,0,0,0,1,1,0,1,1,0,1,1,1,1,0,
    0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,0,
    0,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,
    1,1,1,1,1,1,1,1,0,1,0,1,0,1,1,1,1,1,0,1),
  nrow = 20, byrow = TRUE
)

maze <- raster(maze, xmn = 0.5, xmx = ncol(maze) + 0.5, ymn = 0.5, ymx = nrow(maze) + 0.5)
maze[maze==0] <- NA # 0 makes the formatting cleaner above, but NA is needed for true barriers

# Get info about the shortest path through the maze using gdistance
lcd <- (function() {
  points <- xyFromCell(maze, c(1, 400))

  tr <- transition(maze, function(x) 1/mean(x), 4)
  tr <- geoCorrection(tr)

  list(dist = costDistance(tr, points),
       path = shortestPath(tr, points[1, ], points[2, ], output="SpatialLines"))
})()

# Basic maze layout
maze_plot(maze, "Resistance", vir_col[1])
lines(lcd$path, col = vir_col[2], lw = 3)


## @knitr 1_setup_4
# End of maze
maze_finish <- maze * 0
maze_finish[20, 20] <- 1

maze_plot(maze_finish, "Absorption", vir_col)


## @knitr 1_setup_5
tolerance = sqrt(.Machine$double.eps) # Default tolerance in functions like all.equal()

print(tolerance)


## @knitr 1_setup_6
tr <- list(fun = function(x) 1/mean(x), dir = 4, sym = TRUE)

samc_obj <- samc(maze, maze_finish, tr_args = tr)

start <- locate(samc_obj, data.frame(x = 1, y = 20))
finish <- locate(samc_obj, data.frame(x = 20, y = 1))


## @knitr 1_ttf_1
survive <- survival(samc_obj)

maze_plot(map(samc_obj, survive), "Expected time to finish", viridis(256))


## @knitr 1_ttf_2
survive[start]


## @knitr 1_ttf_3
cond <- cond_passage(samc_obj, dest = finish)

cond[start]


## @knitr 1_pov_1
disp <- dispersal(samc_obj, origin = start)

maze_plot(map(samc_obj, disp), "Probability of Visit", viridis(256))


## @knitr 1_pov_2
# Ideally would use `as.numeric(disp == 1)`, but floating point precision issues force an approximation
disp_sol <- as.numeric(abs(disp - 1) < tolerance)

maze_plot(map(samc_obj, disp_sol), "Solution Using Dispersal()", vir_col)


## @knitr 1_pov_3
disp[start]


## @knitr 1_visit_1
visit <- visitation(samc_obj, origin = start)

maze_plot(map(samc_obj, visit), "Visits Per Cell", viridis(256))


## @knitr 1_visit_2
visit[finish]


## @knitr 1_loc_1
dist <- distribution(samc_obj, origin = start, time = 20)

maze_plot(map(samc_obj, dist), "Location at t=20", col = viridis(256))


## @knitr 1_loc_2
dist <- distribution(samc_obj, origin = start, time = 21)

maze_plot(map(samc_obj, dist), "Location at t=21", viridis(256))


## @knitr 1_occ_1
maze_occ <- maze * 0
maze_occ[1, 1] <- 1

maze_plot(maze_occ, "Occupancy", vir_col)


## @knitr 1_occ_2
# Scenario 1: 3 people start in the maze
maze_occ3 <- maze * 0
maze_occ3[1, 1] <- 3

survival(samc_obj, occ = maze_occ3)


## @knitr 1_occ_3
survival(samc_obj, occ = maze_occ3) / 3



## @knitr 1_occ_4
# Scenario 2: A person starts in each corner of the maze
maze_occ3 <- maze * 0
maze_occ3[1, 1] <- 1
maze_occ3[20, 1] <- 1
maze_occ3[1, 20] <- 1

maze_plot(maze_occ, "Occupancy", vir_col)

survival(samc_obj, occ = maze_occ3)


## @knitr 1_occ_5
survival(samc_obj, occ = maze_occ3) / 3


## @knitr 1_occ_6
dist <- distribution(samc_obj, occ = maze_occ3, time = 17)

# This makes it easier to see how far along the individuals could be
dist <- as.numeric(dist > 0)

maze_plot(map(samc_obj, dist), "Location at t=17", viridis(256))



## @knitr Part2
#
# Part 2 ----
#

## @knitr 2_fid_1
# Intersections determined using a moving window function
maze_ints <- focal(maze, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) > 3}, pad = TRUE)
maze_ints[is.na(maze)] <- NA
maze_ints <- maze_ints * 0.1

maze_plot(maze_ints, "Intersections", vir_col)

## @knitr 2_fid_2
samc_ints <- samc(maze, maze_finish, maze_ints, tr_args = tr)


## @knitr 2_fid_3
# Original results from Part 1
survival(samc_obj)[start]
cond_passage(samc_obj, start, finish)

# Results with fidelity at intersections
survival(samc_ints)[start]
cond_passage(samc_ints, start, finish)


## @knitr 2_fid_4
disp <- dispersal(samc_obj, origin = start)
disp_ints <- dispersal(samc_ints, origin = start)

all.equal(disp, disp_ints)


## @knitr 2_fid_5
visit_orig <- visitation(samc_obj, origin = start)
visit_ints <- visitation(samc_ints, origin = start)

all.equal(visit_orig, visit_ints)

# Let's plot the difference to see if there is a noticeable pattern
visit_diff <- map(samc_obj, visit_ints) - map(samc_obj, visit_orig)
maze_plot(visit_diff, "Visits Per Cell (Difference)", viridis(256))


## @knitr 2_fid_6
# First, let's see which cells changed.
# Ideally would just use `visit_diff > 0`, but floating point precision issues force an approximation
maze_plot(visit_diff > tolerance, "Visits With Non-Zero Difference", vir_col)

# Second, let's see what the percent change is for our non-zero differences.
visit_perc <- (visit_ints - visit_orig) / visit_orig
visit_perc[visit_perc>tolerance]


## @knitr 2_fid_7
dist_ints <- distribution(samc_ints, origin = start, time = 20)
maze_plot(map(samc_ints, dist_ints), "Location at t=20", viridis(256))

dist_ints <- distribution(samc_ints, origin = start, time = 21)
maze_plot(map(samc_ints, dist_ints), "Location at t=21", viridis(256))


## @knitr 2_fid_8
dist_ints <- distribution(samc_ints, origin = start, time = 200)
maze_plot(map(samc_ints, dist_ints), "Location at t=200", viridis(256))

dist_ints <- distribution(samc_ints, origin = start, time = 201)
maze_plot(map(samc_ints, dist_ints), "Location at t=201", viridis(256))


## @knitr 2_fid_9
dist <- distribution(samc_obj, origin = start, time = 200)
maze_plot(map(samc_obj, dist), "Location at t=200", viridis(256))

dist <- distribution(samc_obj, origin = start, time = 201)
maze_plot(map(samc_obj, dist), "Location at t=201", viridis(256))


## @knitr 2_end_1
# Dead ends
maze_ends <- focal(maze, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) == 2}, pad = TRUE)
maze_ends[is.na(maze)] <- NA
maze_ends <- maze_ends * 9 + 1
maze_ends[20, 20] <- 1

maze_plot(maze_ends, "Dead Ends", vir_col)


## @knitr 2_end_2
samc_ends <- samc(maze_ends, maze_finish, tr_args = tr)


## @knitr 2_end_3
# Original results from Part 1
survival(samc_obj)[start]
cond_passage(samc_obj, start, finish)

# Results with dead ends
survival(samc_ends)[start]
cond_passage(samc_ends, start, finish)


## @knitr 2_end_4
disp_ends <- dispersal(samc_ends, origin = start)
maze_plot(map(samc_obj, disp_ends), "Probability of Visit", viridis(256))

visit_ends <- visitation(samc_ends, origin = start)
maze_plot(map(samc_obj, visit_ends), "Visits Per Cell", viridis(256))





## @knitr 2_traps_1
# Traps absorption layer
maze_traps <- maze * 0
maze_traps[17, 3] <- 0.2
maze_traps[1, 9] <- 0.2
maze_traps[6, 20] <- 0.2

maze_plot(maze_traps, "Traps", vir_col)


## @knitr 2_traps_2
abs_total <- maze_finish + maze_traps

samc_traps <- samc(maze, abs_total, tr_args = tr)


## @knitr 2_traps_3
# Original results from Part 1
survival(samc_obj)[start]
cond_passage(samc_obj, start, finish)

# Results with traps
survival(samc_traps)[start]
cond_passage(samc_traps, start, finish)


## @knitr 2_traps_4
survive_traps <- survival(samc_traps)

# Note the updated title from part 1
maze_plot(map(samc_obj, survive_traps), "Expected Time to Absorption", viridis(256))


## @knitr 2_traps_5
disp_traps <- dispersal(samc_traps, origin = start)
maze_plot(map(samc_traps, disp_traps), "Probability of Visit", viridis(256))

visit_traps <- visitation(samc_traps, origin = start)
maze_plot(map(samc_traps, visit_traps), "Visits Per Cell", viridis(256))


## @knitr 2_traps_6
# Ideally, we would just use `as.numeric(disp == 1)`, but we have floating point precision issues here, so we will approximate it
disp_traps_route <- as.numeric(abs(disp_traps - 1) < tolerance)

maze_plot(map(samc_traps, disp_traps_route), "dispersal() == 1", vir_col)




## @knitr 2_add_1
mort_traps <- mortality(samc_traps, origin = start)

maze_plot(map(samc_traps, mort_traps), "Absorption Probability", viridis(256))


## @knitr 2_add_2
mort_traps[mort_traps > 0]

mort_traps[finish]


## @knitr 2_add_3
# Naming the rasters will make things easier and less prone to user error later
names(maze_finish) <- "Finish"
names(maze_traps) <- "Traps"

samc_traps$abs_states <- raster::stack(maze_finish, maze_traps)


## @knitr 2_add_4
mort_traps <- mortality(samc_traps, origin = start)

str(mort_traps)

maze_plot(map(samc_traps, mort_traps$Finish), "Absorption Probability (Finish)", viridis(256))
maze_plot(map(samc_traps, mort_traps$Traps), "Absorption Probability (Traps)", viridis(256))


## @knitr 2_add_5
absorption(samc_traps, origin = start)


## @knitr Part3
#
# Part 3 ----
#
