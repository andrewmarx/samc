# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

# This script is the source of the code for maze example vignettes
# The lines with `@knitr` are for processing by knitr


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

## @knitr maze_ints
# Intersections determined using a moving window function
maze_ints <- focal(maze, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) > 3}, pad = TRUE)
maze_ints[is.na(maze)] <- NA
maze_ints <- maze_ints * 0.1

maze_plot(maze_ints, "Intersections", vir_col)

## @knitr samc_ints
samc_ints <- samc(maze, maze_finish, maze_ints, tr_args = tr)


## @knitr maze_ends
# Dead ends
maze_ends <- focal(maze, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) == 2}, pad = TRUE)
maze_ends[is.na(maze)] <- NA
maze_ends <- maze_ends * 9 + 1
maze_ends[20, 20] <- 1


## @knitr samc_ends
samc_ends <- samc(maze_ends, maze_finish, tr_args = tr)


## @knitr maze_traps
# Traps absorption layer
maze_traps <- maze * 0
maze_traps[17, 3] <- 0.2
maze_traps[1, 9] <- 0.2
maze_traps[6, 20] <- 0.2

## @knitr samc_traps
abs_total <- maze_finish + maze_traps

samc_traps <- samc(maze, abs_total, tr_args = tr)


## @knitr mort_traps
mort_traps <- mortality(samc_traps, origin = start)



## @knitr Part3
#
# Part 3 ----
#
