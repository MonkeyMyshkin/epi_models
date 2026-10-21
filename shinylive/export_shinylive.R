# Export the practical Shiny apps as one static Shinylive site that runs
# entirely in the students' browsers (no R or Shiny server needed). The apps
# share a single copy of the webR runtime, and site/index.html links to them.
#
# Usage, from the Practicals folder:
#   Rscript shinylive/export_shinylive.R
# The site is written to _site/ (one subfolder per app) and must be served over
# http(s); it won't work opened as a file. GitHub Actions builds and publishes it
# on every push to main (.github/workflows/deploy.yml).
# To try it locally: Rscript -e 'httpuv::runStaticServer("_site")'
#
# The apps avoid deSolve (no working webR build) by solving the ODEs with a
# base-R rk4() function, so every package they use is available in webR.

apps    <- c("SIRmodel", "SEIRmodel", "SISmodel")
out_dir <- "_site"

unlink(out_dir, recursive = TRUE)
for(app in apps)
{
    shinylive::export(app, out_dir, subdir = app)
}
file.copy("site/index.html", file.path(out_dir, "index.html"))
