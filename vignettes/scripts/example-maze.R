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
plot_maze <- function(map, title, colors) {
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
maze_res = matrix(
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

maze_res <- raster(maze_res, xmn = 0.5, xmx = ncol(maze_res) + 0.5, ymn = 0.5, ymx = nrow(maze_res) + 0.5)
maze_res[maze_res==0] <- NA # 0 makes the formatting cleaner above, but NA is needed for true barriers

# Get info about the shortest path through the maze using gdistance
lcd <- (function() {
  points <- xyFromCell(maze_res, c(1, 400))

  tr <- transition(maze_res, function(x) 1/mean(x), 4)
  tr <- geoCorrection(tr)

  list(dist = costDistance(tr, points),
       path = shortestPath(tr, points[1, ], points[2, ], output="SpatialLines"))
})()

# Basic maze layout
plot_maze(maze_res, "Resistance", vir_col[1])
lines(lcd$path, col = vir_col[2], lw = 3)


## @knitr 1_setup_4
# End of maze
maze_finish <- maze_res * 0
maze_finish[20, 20] <- 1

plot_maze(maze_finish, "Absorption", vir_col)


## @knitr 1_setup_5
tolerance = sqrt(.Machine$double.eps) # Default tolerance in functions like all.equal()

print(tolerance)


## @knitr 1_setup_6
tr <- list(fun = function(x) 1/mean(x), dir = 4, sym = TRUE)

maze_samc <- samc(maze_res, maze_finish, tr_args = tr)

maze_origin <- locate(maze_samc, data.frame(x = 1, y = 20))
maze_dest <- locate(maze_samc, data.frame(x = 20, y = 1))


## @knitr 1_ttf_1
maze_surv <- survival(maze_samc)

plot_maze(map(maze_samc, maze_surv), "Expected time to finish", viridis(256))


## @knitr 1_ttf_2
maze_surv[maze_origin]


## @knitr 1_ttf_3
maze_cond <- cond_passage(maze_samc, dest = maze_dest)

maze_cond[maze_origin]


## @knitr 1_pov_1
maze_disp <- dispersal(maze_samc, origin = maze_origin)

plot_maze(map(maze_samc, maze_disp), "Probability of Visit", viridis(256))


## @knitr 1_pov_2
# Ideally would use `as.numeric(maze_disp == 1)`, but floating point precision issues force an approximation
maze_disp_sol <- as.numeric(abs(maze_disp - 1) < tolerance)

plot_maze(map(maze_samc, maze_disp_sol), "Solution Using Dispersal()", vir_col)


## @knitr 1_pov_3
maze_disp[maze_origin]


## @knitr 1_visit_1
maze_visit <- visitation(maze_samc, origin = maze_origin)

plot_maze(map(maze_samc, maze_visit), "Visits Per Cell", viridis(256))


## @knitr 1_visit_2
maze_visit[maze_dest]


## @knitr 1_loc_1
maze_dist <- distribution(maze_samc, origin = maze_origin, time = 20)

plot_maze(map(maze_samc, maze_dist), "Location at t=20", col = viridis(256))


## @knitr 1_loc_2
maze_dist <- distribution(maze_samc, origin = maze_origin, time = 21)

plot_maze(map(maze_samc, maze_dist), "Location at t=21", viridis(256))


## @knitr 1_occ_1
maze_occ <- maze_res * 0
maze_occ[1, 1] <- 1

plot_maze(maze_occ, "Occupancy", vir_col)


## @knitr 1_occ_2
survival(maze_samc, occ = maze_occ)

maze_surv[maze_origin]


## @knitr 1_occ_3
# Scenario 1: 3 people start in the maze
maze_occ3 <- maze_res * 0
maze_occ3[1, 1] <- 3

survival(maze_samc, occ = maze_occ3)


## @knitr 1_occ_4
survival(maze_samc, occ = maze_occ3) / 3


## @knitr 1_occ_5
# Scenario 2: A person starts in each corner of the maze
maze_occ3 <- maze_res * 0
maze_occ3[1, 1] <- 1
maze_occ3[20, 1] <- 1
maze_occ3[1, 20] <- 1

plot_maze(maze_occ, "Occupancy", vir_col)

survival(maze_samc, occ = maze_occ3)


## @knitr 1_occ_6
survival(maze_samc, occ = maze_occ3) / 3


## @knitr 1_occ_7
maze_occ3_dist <- distribution(maze_samc, occ = maze_occ3, time = 17)

# This makes it easier to see how far along the individuals could be
maze_occ3_dist <- as.numeric(maze_occ3_dist > 0)

plot_maze(map(maze_samc, maze_occ3_dist), "Location at t=17", viridis(256))



## @knitr Part2
#
# Part 2 ----
#

## @knitr 2_fid_1
# Intersections determined using a moving window function
ints_res <- focal(maze_res, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) > 3}, pad = TRUE)
ints_res[is.na(maze_res)] <- NA
ints_res <- ints_res * 0.1

plot_maze(ints_res, "Intersections", vir_col)

## @knitr 2_fid_2
ints_samc <- samc(maze_res, maze_finish, ints_res, tr_args = tr)


## @knitr 2_fid_3
# Original results from Part 1
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, maze_origin, maze_dest)

# Results with fidelity at intersections
survival(ints_samc)[maze_origin]
cond_passage(ints_samc, maze_origin, maze_dest)


## @knitr 2_fid_4
ints_disp <- dispersal(ints_samc, origin = maze_origin)

all.equal(maze_disp, ints_disp)


## @knitr 2_fid_5
ints_visit <- visitation(ints_samc, origin = maze_origin)

all.equal(maze_visit, ints_visit)

# Let's plot the difference to see if there is a noticeable pattern
visit_diff <- map(maze_samc, ints_visit) - map(maze_samc, maze_visit)
plot_maze(visit_diff, "Visits Per Cell (Difference)", viridis(256))


## @knitr 2_fid_6
# First, let's see which cells changed.
# Ideally would just use `visit_diff > 0`, but floating point precision issues force an approximation
plot_maze(visit_diff > tolerance, "Visits With Non-Zero Difference", vir_col)

# Second, let's see what the percent change is for our non-zero differences.
visit_perc <- (ints_visit - maze_visit) / maze_visit
visit_perc[visit_perc>tolerance]


## @knitr 2_fid_7
ints_dist <- distribution(ints_samc, origin = maze_origin, time = 20)
plot_maze(map(ints_samc, ints_dist), "Location at t=20", viridis(256))

ints_dist <- distribution(ints_samc, origin = maze_origin, time = 21)
plot_maze(map(ints_samc, ints_dist), "Location at t=21", viridis(256))


## @knitr 2_fid_8
ints_dist <- distribution(ints_samc, origin = maze_origin, time = 200)
plot_maze(map(ints_samc, ints_dist), "Location at t=200", viridis(256))

ints_dist <- distribution(ints_samc, origin = maze_origin, time = 201)
plot_maze(map(ints_samc, ints_dist), "Location at t=201", viridis(256))


## @knitr 2_fid_9
maze_dist <- distribution(maze_samc, origin = maze_origin, time = 200)
plot_maze(map(maze_samc, maze_dist), "Location at t=200", viridis(256))

maze_dist <- distribution(maze_samc, origin = maze_origin, time = 201)
plot_maze(map(maze_samc, maze_dist), "Location at t=201", viridis(256))


## @knitr 2_end_1
# Dead ends
ends_res <- focal(maze_res, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) == 2}, pad = TRUE)
ends_res[is.na(maze_res)] <- NA
ends_res <- ends_res * 9 + 1
ends_res[20, 20] <- 1

plot_maze(ends_res, "Dead Ends", vir_col)


## @knitr 2_end_2
ends_samc <- samc(ends_res, maze_finish, tr_args = tr)


## @knitr 2_end_3
# Original results from Part 1
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, maze_origin, maze_dest)

# Results with dead ends
survival(ends_samc)[maze_origin]
cond_passage(ends_samc, maze_origin, maze_dest)


## @knitr 2_end_4
ends_disp <- dispersal(ends_samc, origin = maze_origin)
plot_maze(map(maze_samc, ends_disp), "Probability of Visit", viridis(256))

ends_visit <- visitation(ends_samc, origin = maze_origin)
plot_maze(map(maze_samc, ends_visit), "Visits Per Cell", viridis(256))





## @knitr 2_traps_1
# Traps absorption layer
maze_traps <- maze_res * 0
maze_traps[17, 3] <- 0.2
maze_traps[1, 9] <- 0.2
maze_traps[6, 20] <- 0.2

plot_maze(maze_traps, "Traps", vir_col)


## @knitr 2_traps_2
maze_abs_total <- maze_finish + maze_traps

traps_samc <- samc(maze_res, maze_abs_total, tr_args = tr)


## @knitr 2_traps_3
# Original results from Part 1
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, maze_origin, maze_dest)

# Results with traps
survival(traps_samc)[maze_origin]
cond_passage(traps_samc, maze_origin, maze_dest)


## @knitr 2_traps_4
traps_surv <- survival(traps_samc)

# Note the updated title from part 1
plot_maze(map(maze_samc, traps_surv), "Expected Time to Absorption", viridis(256))


## @knitr 2_traps_5
traps_disp <- dispersal(traps_samc, origin = maze_origin)
plot_maze(map(traps_samc, traps_disp), "Probability of Visit", viridis(256))

traps_visit <- visitation(traps_samc, origin = maze_origin)
plot_maze(map(traps_samc, traps_visit), "Visits Per Cell", viridis(256))


## @knitr 2_traps_6
# Ideally, we would just use `as.numeric(traps_disp == 1)`, but we have floating point precision issues here, so we will approximate it
traps_disp_route <- as.numeric(abs(traps_disp - 1) < tolerance)

plot_maze(map(traps_samc, traps_disp_route), "dispersal() == 1", vir_col)


## @knitr 2_add_1
traps_mort <- mortality(traps_samc, origin = maze_origin)

plot_maze(map(traps_samc, traps_mort), "Absorption Probability", viridis(256))


## @knitr 2_add_2
traps_mort[traps_mort > 0]

traps_mort[maze_dest]


## @knitr 2_add_3
# Naming the rasters will make things easier and less prone to user error later
names(maze_finish) <- "Finish"
names(maze_traps) <- "Traps"

traps_samc$abs_states <- raster::stack(maze_finish, maze_traps)


## @knitr 2_add_4
traps_mort_dec <- mortality(traps_samc, origin = maze_origin)

str(traps_mort_dec)

plot_maze(map(traps_samc, traps_mort_dec$Finish), "Absorption Probability (Finish)", viridis(256))
plot_maze(map(traps_samc, traps_mort_dec$Traps), "Absorption Probability (Traps)", viridis(256))


## @knitr 2_add_5
absorption(traps_samc, origin = maze_origin)


## @knitr Part3
#
# Part 3 ----
#

## @knitr 3_short_1
# Create a copy and add a shortcut
short_res <- maze_res
short_res[16, 6] <- 10

# Get info about the shortest path through the new maze using gdistance
lcd2 <- (function() {
  points <- xyFromCell(short_res, c(1, 400))

  tr <- transition(short_res, function(x) 1/mean(x), 4)
  tr <- geoCorrection(tr)

  list(dist = costDistance(tr, points),
       path = shortestPath(tr, points[1, ], points[2, ], output="SpatialLines"))
})()

plot_maze(short_res, "Shortcut Maze", vir_col)
lines(lcd2$path, col = vir_col[2], lw = 3)


## @knitr 3_short_2
# Let's see what the difference in distance is
lcd2$dist - lcd$dist


## @knitr 3_short_3
# Our old absorption layer does not quite match our new resistance layer, so make a new one
short_finish <- short_res * 0
short_finish[20, 20] <- 1


## @knitr 3_short_4
short_samc <- samc(short_res, short_finish, tr_args = tr)

# Important: we have to rerun locate()
short_origin <- locate(short_samc, data.frame(x = 1, y = 20))
short_dest <- locate(short_samc, data.frame(x = 20, y = 1))


## @knitr 3_short_5
short_surv <- survival(short_samc)

plot_maze(map(short_samc, short_surv), "Expected time to finish (Shortcut Maze)", viridis(256))


## @knitr 3_short_6
# Expected time to finish from the start
short_surv[maze_origin]

# The difference from our original maze
short_surv[maze_origin] - maze_surv[maze_origin]


## @knitr 3_short_7
short_cond <- cond_passage(short_samc, dest = short_dest)
short_cond[maze_origin]


## @knitr 3_short_8
short_disp <- dispersal(short_samc, origin = short_origin)

plot_maze(map(short_samc, short_disp), "Probability of Visit (Shortcut Maze)", viridis(256))


## @knitr 3_short_9
# Ideally, we would just use `as.numeric(short_disp == 1)`, but we have floating point precision issues here, so we will approximate it
short_disp_sol <- as.numeric(abs(short_disp - 1) < tolerance)

plot_maze(map(short_samc, short_disp_sol), "Partial solution (Shortcut Maze)", vir_col)


## @knitr 3_combine_1
# Combine our previous resistance layers
all_res <- max(stack(short_res, ints_res, ends_res), na.rm = TRUE)

# For absorption, all we need is an updated version of our traps raster
all_traps <- maze_traps
all_traps[16, 6] <- 0

# Total absorption
all_abs_total <- short_finish + all_traps


# If we had more variety in our resistance values we would want more colors
plot_maze(all_res, "Final Maze", vir_col)

# Plot the traps raster
plot_maze(all_traps, "Final Maze Traps", vir_col)


## @knitr 3_combine_2
all_samc <- samc(all_res, all_abs_total, tr_args = tr)

# We can actually reuse the short_res locations in this case, but let's make new ones anyway
all_start <- locate(all_samc, data.frame(x = 1, y = 20))
all_finish <- locate(all_samc, data.frame(x = 20, y = 1))


## @knitr 3_combine_3
all_surv <- survival(all_samc)

# Note the updated title from part 1
plot_maze(map(all_samc, all_surv), "Expected Time to Absorption", viridis(256))


## @knitr 3_combine_4
# Original results (Part 1)
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, maze_origin, maze_dest)

# Results with traps (Part 2)
survival(traps_samc)[maze_origin]
cond_passage(traps_samc, maze_origin, maze_dest)

# Results with a shortcut
survival(short_samc)[short_origin]
cond_passage(short_samc, short_origin, short_dest)

# Results with everything
survival(all_samc)[all_start]
cond_passage(all_samc, all_start, all_finish)


## @knitr 3_combine_5
traps_mort[traps_mort > 0]

all_mort <- mortality(all_samc, origin = all_start)
all_mort[all_mort > 0]
