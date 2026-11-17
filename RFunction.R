library(moveapps)
library(ctmm)
library(sf)

rFunction <- function(data, isopleth_levels) {
  lvls <- parse_levels(isopleth_levels)

  models <- data[[1]]
  telemetry <- data[[2]]

  logger.info(
    paste0(
      "Estimating aKDE home ranges for ", length(telemetry),
      " animal(s) at isopleth level(s) ", paste(lvls, collapse = ", ")
    )
  )
  
  hr <- akde(telemetry, models)
  
  # Artefact: home range contours of all animals at all levels
  contours <- purrr::list_rbind(
    purrr::imap(hr, function(rng, id) animal_contours(rng, id, lvls))
  )
  
  sf::st_write(
    contours,
    appArtifactPath("homerange.gpkg"),
    append = FALSE,
    quiet = TRUE
  )

  # Artefact: home range areas
  areas <- purrr::list_rbind(
    purrr::imap(hr, function(rng, id) animal_summary(rng, id, lvls))
  )
  
  write.csv(areas, appArtifactPath("akde_summary.csv"), row.names = FALSE)

  write_pdf(appArtifactPath("akde_plots.pdf"), telemetry, hr, lvls)

  # Artefact: one UD raster per animal. These do not depend on the isopleth level.
  targetDirUDs <- tempfile("uds")
  dir.create(targetDirUDs)
  
  purrr::walk(
    names(hr),
    function(id) {
      writeRaster(
        hr[[id]],
        file.path(targetDirUDs, paste0(id, ".tif")),
        overwrite = TRUE
      )
    }
  )
  
  zip::zip(
    zipfile = appArtifactPath("uds.zip"),
    files = list.files(targetDirUDs, full.names = TRUE, pattern = "tif$"),
    mode = "cherry-pick"
  )

  # List elements need names for cargo agent analyzer. Note that some apps
  # may consume this output with positional matching though.
  list(model = models, uds = hr, data = telemetry)
}


# Convert isopleth level setting (comma delimited string) into a sorted
# vector of unique numbers. Stops with a message if any value is not a number
# strictly between 0 and 1.
parse_levels <- function(isopleth_levels) {
  if (is.null(isopleth_levels) || all(trimws(isopleth_levels) == "")) {
    stop("No isopleth level given. Enter one or more levels between 0 and 1, ",
         "separated by commas (e.g. \"0.5, 0.95\").", call. = FALSE)
  }
  
  parts <- trimws(unlist(strsplit(as.character(isopleth_levels), ",")))
  levels <- suppressWarnings(as.numeric(parts))
  
  bad <- parts[is.na(levels) | levels <= 0 | levels >= 1]
  
  if (length(bad) > 0) {
    stop("Invalid isopleth level(s): ", paste0("\"", bad, "\"", collapse = ", "),
         ". Each level must be a number between 0 and 1, separated by commas ",
         "(e.g. \"0.5, 0.95\"). For a 95% home range, enter 0.95, not 95.",
         call. = FALSE)
  }
  
  sort(unique(levels))
}

# Generate contour polygons for the home range at the given isopleth level.
# This ensures that internal holes within the home range contours are
# correctly left unfilled. At the time of writing, there is a bug in ctmm
# (1.3.x) that does not mark these holes correctly, meaning that the entire
# region within the external contour is filled.
ud_contours <- function(ud, level.UD) {
  spdf <- ctmm::SpatialPolygonsDataFrame.UD(ud, level.UD = level.UD)

  # Use even-odd rule to identify holes. Contour rings are nested by definition,
  # so the symmetric difference of all rings keeps the area inside an odd
  # number of rings and drops the areas inside an even number (holes).
  geometries <- lapply(spdf@polygons, function(contour) {
    rings <- lapply(contour@Polygons, function(ring) sf::st_polygon(list(ring@coords)))
    geom <- Reduce(sf::st_sym_difference, rings)

    # Rings that touch can leave zero-width line fragments in the result.
    # Keep only the polygon parts.
    if (inherits(geom, "GEOMETRYCOLLECTION")) {
      geom <- sf::st_union(sf::st_collection_extract(sf::st_sfc(geom), "POLYGON"))[[1]]
    }

    sf::st_cast(geom, "MULTIPOLYGON")
  })

  # Rebuild output to mirror `ctmm::SpatialPolygonsDataFrame.UD()` structure.
  # Contour names are the ctmm polygon `ID` values.
  sf::st_sf(
    name = vapply(spdf@polygons, function(contour) contour@ID, character(1)),
    geometry = sf::st_sfc(geometries, crs = sf::st_crs(spdf@proj4string))
  )
}

# Get contours of one animal's UD at each isopleth level. Each row represents
# a contour for an animal, level, and estimate type (estimate or CI)
animal_contours <- function(ud, id, levels) {
  lvl_contours <- purrr::map(
    levels,
    function(lvl) {
      contours <- ud_contours(ud, level.UD = lvl)
      
      contours$id <- id
      contours$level <- lvl
      contours$estimate <- sub(".* ", "", contours$name)
      contours[, c("id", "level", "estimate", "name")]
    }
  )
  
  purrr::list_rbind(lvl_contours)
}

# Get home range area of one animal at each isopleth level, in square meters
# so that values are comparable across animals and levels.
animal_summary <- function(ud, id, levels) {
  lvl_summaries <- purrr::map(
    levels,
    function(lvl) {
      ci <- summary(ud, level.UD = lvl, units = FALSE)$CI
      
      data.frame(
        id = id,
        level = lvl,
        parameter = rownames(ci),
        low = ci[, "low"],
        est = ci[, "est"],
        high = ci[, "high"],
        row.names = NULL
      )
    }
  )
  
  purrr::list_rbind(lvl_summaries)
}

# Colors for the contours at each isopleth level. Viridis palette starting from
# teal to dark purple. Stops at teal because lighter viridis palette colors
# are hard to make out on white background.
contour_colors <- function(n) {
  grDevices::hcl.colors(100, "viridis")[round(seq(55, 5, length.out = n))]
}

# Legend shown on each animal's page, with colors for each isopleth contour
# and the area/95% CI for each isopleth range.
area_legend <- function(ud, levels) {
  ci <- lapply(levels, function(level) summary(ud, level.UD = level)$CI)

  unit <- sub(".*\\((.*)\\)", "\\1", rownames(ci[[1]])[1])
  fmt <- function(x) format(signif(x, 3), big.mark = ",")
  cols <- contour_colors(length(levels))

  par(mar = c(0, 0, 4.1, 0))
  plot.new()
  cex <- 0.9
  line_height <- 1.5 * strheight("M", cex = cex)
  text_x <- 0.13

  y <- 1
  text(0, y, "Home range area", adj = c(0, 1), cex = cex, font = 2)
  y <- y - line_height
  text(0, y, paste0("(", unit, " with 95% CI)"), adj = c(0, 1), cex = cex)

  for (i in seq_along(levels)) {
    x <- ci[[i]]
    y <- y - 2 * line_height
    segments(0, y - line_height / 3, 0.1, y - line_height / 3, col = cols[i], lwd = 3)
    text(text_x, y, paste0(levels[i] * 100, "% range"), adj = c(0, 1), cex = cex)
    y <- y - line_height
    text(text_x, y, paste0(fmt(x[1, "est"]), " [", fmt(x[1, "low"]), " - ", fmt(x[1, "high"]), "]"),
         adj = c(0, 1), cex = cex)
  }
}

# Legend of animal colours for the side panel of an overview page
animal_legend <- function(ids, cols) {
  par(mar = c(0, 0, 4.1, 0))
  plot.new()
  legend("topleft", legend = ids, col = cols, pch = 19, bty = "n",
         ncol = ceiling(length(ids) / 25),
         cex = if (length(ids) > 25) 0.7 else 1)
}

# PDF summary artifact. One page per isopleth with all animals together,
# then individual animal pages with all isopleths overlaid onto the single
# home range density gradient.
write_pdf <- function(file, telemetry, hr, levels) {
  animal_cols <- grDevices::hcl.colors(length(hr), palette = "Dark 3")
  plot_mar <- c(5.1, 4.1, 4.1, 2.1)

  pdf(file, width = 11, height = 7.5)
  on.exit(dev.off())
  layout(matrix(1:2, nrow = 1), widths = c(3, 1))

  for (level in levels) {
    par(mar = plot_mar)
    plot(telemetry, UD = hr, level.UD = level, col = animal_cols,
         col.UD = animal_cols, col.level = animal_cols, col.grid = NA)
    title(paste0("All animals: ", level * 100, "% home range"))
    animal_legend(names(hr), animal_cols)
  }

  for (id in names(hr)) {
    par(mar = plot_mar)
    # The locations are drawn in a second call. ctmm sizes the location error
    # circles from level.UD, and with several levels it mixes their values up
    # across the locations (with R warnings). A second call with the default
    # single level draws the circles as in any standard ctmm plot.
    plot(telemetry[[id]], UD = hr[[id]], level.UD = levels, col = NA, error = FALSE,
         col.UD = "grey45", col.level = contour_colors(length(levels)), lwd.level = 2)
    plot(telemetry[[id]], col = "black", add = TRUE)
    title(id)
    area_legend(hr[[id]], levels)
  }
}
