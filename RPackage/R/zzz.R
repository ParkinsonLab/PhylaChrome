PhylaChromeMain_mod <- NULL

.onLoad <- function(libname, pkgname) {
  reticulate::py_require(c("attrs", "chromato", "pandas", "multitax", "scipy"))

  # Locate the Python script inside the installed package environment
  script_path <- system.file("python", "PhylaChromeMain.py", package = pkgname)

  # Fallback helper for local interactive testing (devtools::load_all)
  if (script_path == "") {
    script_path <- file.path(getwd(), "inst", "python", "PhylaChromeMain.py")
  }

  # Delay module loading until Python is initialized by the user/system
  PhylaChromeMain_mod <<- reticulate::import_from_path(
    module = "PhylaChromeMain",
    path = dirname(script_path),
    convert = TRUE, # Automatically converts data types between R and Python
    delay_load = TRUE
  )
}
