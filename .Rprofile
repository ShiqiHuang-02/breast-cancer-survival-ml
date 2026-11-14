# Startup: make the project-local package library available (lightweight renv substitute).
local({
  lib <- file.path(getwd(), "lib")
  if (dir.exists(lib)) {
    .libPaths(c(lib, .libPaths()))
  }
})
