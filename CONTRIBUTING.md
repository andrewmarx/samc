# Contributing

### Package Development

Useful references for creating packages:
> http://r-pkgs.had.co.nz/
>
> https://support.rstudio.com/hc/en-us/sections/200130627-Package-Development
>
> http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html
>
> https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/`


### Code Style

All code in the package should conform to Hadley Wickham's tidyverse style guide. Please do not commit changes to the repository until your code conforms to it.

https://style.tidyverse.org/

http://adv-r.had.co.nz/Style.html (A more brief version)

Note that major code sections in the example are delimited with a multiline comment. The first text element in this comment is a brief title followed by four trailing dashes `----`, which Rstudio uses to generate its code outline

Several settings in Rstudio can be enabled to help with ensuring consistent styling.

> *Tools -> Global Options -> Code -> Saving*
>
> Check the "Ensure the source files end with newline" and "Strip trailing horizontal whitespace when saving" options.

> *Tools -> Global Options -> Code -> Diagnostics*
>
> All of the options under **R Diagnostics** are useful. The last one will help in spotting style issues.

### Versioning

Package releases should conform to semantic versioning rules:

https://semver.org/
