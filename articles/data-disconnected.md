# Disconnected Data

*Note: this vignette is far from complete and only contains a brief
discussion (originally from an email) of the topic in the context of
rasters. It will be rewritten and expanded in future releases.*

## Rasters

Imagine you have a raster with 2 islands, and that water separating the
islands in the raster is represented with `NA`’s. These two islands
would be disconnected and will result in
`"Warning: Input contains disconnected regions. This does not work with the cond_passage() metric."`.
The warning isn’t an issue if you don’t need the
[`cond_passage()`](https://andrewmarx.github.io/samc/reference/cond_passage.md)
function.

Now, each island (aka region) must have at least one pixel with a
non-zero absorbing value. The math algorithms don’t work without that.
If a region does not have at least one non-zero absorbing value, then
the [`samc()`](https://andrewmarx.github.io/samc/reference/samc.md)
function will result in
`Error: All disconnected regions must have at least one non-zero absorption value`

Now the tricky part; an island/region can be as small as a single pixel.
So, if you’re not careful with your raster, you might have very small,
isolated patches as small as a pixel floating around somewhere. If
you’re using just 4 directions for the transition function, they could
even be immediately diagonal to another patch.

The easiest way to spot these is to pick one of your rasters and set all
your non-`NA` values to the same value, like 1. Then turn all the `NA`’s
to a different value, like 0. Then plot it. Code:

``` r
raster[!is.na(raster[])) <- 1 
raster[is.na(raster[])) <- 0
plot(raster)
```

If your raster is particularly large, then you might need to save it to
a high-res image to spot solo pixels in an external picture editor or
viewer.
