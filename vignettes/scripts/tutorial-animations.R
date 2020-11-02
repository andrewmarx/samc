# This script is the source of the code for the animations vignette and is also
# used to generate the output figures. It only needs to be run directly if changes
# to the figures are desired

## @knitr setup
# First step is to load the libraries. Not all of these libraries are stricly
# needed; some are used for convenience and visualization for this tutorial.
library("samc")
library("raster")
library("ggplot2")
library("viridis")
library("gifski")
library("gganimate")


# "Load" the data. In this case we are using data built into the package.
# In practice, users will likely load raster data using the raster() function
# from the raster package.
res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data


# Create a samc object using the resistance and absorption data. We use the
# recipricol of the arithmetic mean for calculating the transition matrix. Note,
# the input data here are matrices, not RasterLayers. If using RasterLayers, the
# `latlon` parameter must be set.
samc_obj <- samc(res_data, abs_data, tr_fun = function(x) 1/mean(x))


# Calculate the probabilities of where an individual starting at specific
# location will be for varying time steps. The starting location is going to
# be cell 1 in the landscape, which is the first non-NA cell going in a
# left-to-right then top-to-bottom order.
time_steps <- ((1:50)*2) ^ 2

dist_list <- distribution(samc_obj, origin = 1, time = time_steps)
dist_map <- map(samc_obj, dist_list)


## @knitr gifski
png_path <- file.path(tempdir(), "frame%03d.png")
png(png_path, width = 6, height = 3, units = "in", res = 100)
par(ask = FALSE)

for (ts in time_steps) {
  name <- as.character(ts)
  plot(dist_map[[name]], main = paste("Individual Location at time step ", name), xlab = "x", ylab = "y", col = viridis(256))
}

dev.off()
png_files <- sprintf(png_path, 1:length(time_steps))
gif_file <- tempfile(fileext = ".gif")
gifski(png_files, gif_file, delay = 0.1, progress = FALSE)
unlink(png_files)
utils::browseURL(gif_file)

## @knitr gifski-save
file.copy(gif_file, "vignettes/img/gifski.gif", overwrite = TRUE)

## @knitr gganimate
# Create an empty dataframe to hold all the data from all the plots
dist_df <- data.frame(x = numeric(0), y = numeric(0), layer = numeric(0), steps = numeric(0))

for (ts in time_steps) {
  name <- as.character(ts)
  dist <- as.data.frame(dist_map[[name]], xy = TRUE, na.rm = TRUE)
  dist$steps <- ts

  dist_df <- rbind(dist_df, dist)
}


# Create the animation. Unfortunately, there does not appear to be a way to
# adjust the color scale dynamically across frames using gganimate at this time
anim <- ggplot(dist_df, aes(x = x, y = y)) +
  geom_raster(aes(fill = layer)) +
  transition_manual(steps) +
  scale_fill_viridis(limits = c(0, max(dist_df$layer))) +
  ggtitle("Individual Location at {current_frame}") +
  coord_equal() +
  theme_bw()

animate(anim, duration = 5, height = 2, width = 6, units = "in", res = 150)

## @knitr gganimate-save
anim_save("vignettes/img/gganimate.gif")
