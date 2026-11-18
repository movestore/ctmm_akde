library(moveapps)
library(ctmm)
library(sf)
library(leaflet)

rFunction <- function(data, isopleth_levels) {
  lvls <- parse_levels(isopleth_levels)

  models <- data[[1]]
  telemetry <- data[[2]]

  # Time since `start`, for the log, e.g. "3.4 secs" or "2.1 mins"
  start <- Sys.time()
  elapsed <- function() format(round(Sys.time() - start, 1))

  logger.info(
    paste0(
      "Estimating aKDE home ranges for ", length(telemetry),
      " animal(s) at isopleth level(s) ", paste(lvls, collapse = ", "),
      ". This can take a while for many animals or long tracks."
    )
  )
  
  hr <- akde(telemetry, models)
  logger.info(paste0("Home ranges estimated after ", elapsed()))

  # akde() returns the model instead of a UD for animals whose model is not
  # range resident. Remove these from the output and warn.
  has_range <- purrr::map_lgl(hr, function(x) inherits(x, "UD"))
  
  if (!any(has_range)) {
    stop("Could not estimate home range for any of the provided animals. ",
         "aKDE is not defined for movement models that are not range-resident ",
         "or tracks that are too short to cover the range.",
         call. = FALSE)
  }

  if (!all(has_range)) {
    logger.warn(
      paste0(
        "Could not estimate home range for animal(s): ",
        paste(names(hr)[!has_range], collapse = ", "), ".\n",
        "aKDE is not defined for movement models that are not range-resident ",
        "or tracks that are too short to cover the range.\n",
        "These animals have been removed from the output."
      )
    )
    hr <- hr[has_range]
    models <- models[has_range]
    telemetry <- telemetry[has_range]
  }

  # Artefact: home range contours of all animals at all levels
  logger.info("Writing home range contours (homerange.gpkg)")
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
  logger.info("Writing home range areas (akde_summary.csv)")
  areas <- purrr::list_rbind(
    purrr::imap(hr, function(rng, id) animal_summary(rng, id, lvls))
  )
  
  write.csv(areas, appArtifactPath("akde_summary.csv"), row.names = FALSE)

  logger.info("Writing plots (akde_plots.pdf)")
  write_pdf(appArtifactPath("akde_plots.pdf"), telemetry, hr, lvls)

  logger.info("Writing interactive map (akde_map.html)")
  write_map(appArtifactPath("akde_map.html"), contours, hr, lvls)

  # Artefact: one UD raster per animal. These do not depend on the isopleth level.
  logger.info("Writing UD rasters (uds.zip)")
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
  logger.info(paste0("All artefacts written after ", elapsed()))

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

# Colour for each animal, used in the PDF overview pages and on the map
animal_colors <- function(n) {
  grDevices::hcl.colors(n, palette = "Dark 3")
}

# Legend shown on each animal's page, with colors for each isopleth contour
# and the area/95% CI for each isopleth range.
area_legend <- function(ud, levels) {
  ci <- lapply(levels, function(level) summary(ud, level.UD = level)$CI)

  fmt <- function(x) format(signif(x, 3), big.mark = ",")
  cols <- contour_colors(length(levels))

  units <- vapply(ci, function(x) sub(".*\\((.*)\\)", "\\1", rownames(x)[1]), character(1))
  range_labels <- paste0(levels * 100, "% range (", units, ")")

  par(mar = c(0, 0, 4.1, 0))
  plot.new()
  text_x <- 0.13
  # Shrink the text if the widest label (e.g. "99.5% range (square
  # kilometers)") would not fit the panel
  cex <- min(0.9, 0.9 * 0.95 * (1 - text_x) / max(strwidth(range_labels, cex = 0.9)))
  line_height <- 1.5 * strheight("M", cex = cex)

  y <- 1
  text(0, y, "Home range area", adj = c(0, 1), cex = cex, font = 2)
  y <- y - line_height
  text(0, y, "(with 95% CI)", adj = c(0, 1), cex = cex)

  for (i in seq_along(levels)) {
    x <- ci[[i]]
    y <- y - 2 * line_height
    segments(0, y - line_height / 3, 0.1, y - line_height / 3, col = cols[i], lwd = 3)
    text(text_x, y, range_labels[i], adj = c(0, 1), cex = cex)
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
  animal_cols <- animal_colors(length(hr))
  plot_mar <- c(5.1, 4.1, 4.1, 2.1)

  pdf(file, width = 11, height = 7.5)
  on.exit(dev.off())
  layout(matrix(1:2, nrow = 1), widths = c(3, 1))

  for (level in levels) {
    par(mar = plot_mar)
    plot(telemetry, UD = hr, level.UD = level, col = animal_cols,
         col.UD = grDevices::adjustcolor(animal_cols, alpha.f = 0.8),
         col.level = animal_cols,
         col.grid = NA)
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

# Interactive map as standalone HTML file. Shows all animals at all requested
# isopleth levels. Animals can be toggled on and off.
write_map <- function(file, contours, hr, levels) {
  ids <- names(hr)
  cols <- stats::setNames(animal_colors(length(ids)), ids)
  contours <- sf::st_transform(contours, 4326)

  # Round the coordinates to 5 decimals (about 1 m), which is plenty for the
  # map and roughly halves the size of the HTML file. Points that become
  # identical to the previous point are dropped.
  round_ring <- function(ring) {
    ring <- round(ring, digits = 5)
    ring[c(TRUE, rowSums(abs(diff(ring))) > 0), , drop = FALSE]
  }
  
  sf::st_geometry(contours) <- sf::st_sfc(
    lapply(
      sf::st_geometry(contours), 
      function(geom) {
        sf::st_multipolygon(
          lapply(geom, function(polygon) lapply(polygon, round_ring))
        )
      }
    ),
    crs = 4326
  )

  # Hover labels for each animal's region
  contours$label <- paste0(contours$id, ", ", contours$level * 100, "% range")

  # Starting map view
  bbox <- unname(
    sf::st_bbox(
      contours[contours$estimate == "est" & contours$level == max(levels), ]
    )
  )

  m <- leaflet(options = leafletOptions(zoomSnap = 0.25, zoomDelta = 0.5)) |>
    addProviderTiles(providers$Esri.WorldTopoMap, group = "World Topo Map") |>
    addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") |>
    fitBounds(bbox[1], bbox[2], bbox[3], bbox[4])

  # Build layer group for each animal, so animals can be toggled on/off
  for (level in sort(levels, decreasing = TRUE)) {
    for (id in ids) {
      this <- contours[contours$id == id & contours$level == level, ]
      m <- m |>
        addPolygons(
          data = this[this$estimate != "est", ],
          group = id,
          label = ~label,
          color = cols[[id]],
          weight = 1,
          opacity = 0.8,
          dashArray = "4",
          fill = FALSE
        ) |>
        addPolygons(
          data = this[this$estimate == "est", ],
          group = id,
          label = ~label,
          color = cols[[id]],
          weight = 2,
          opacity = 1,
          fillColor = cols[[id]],
          fillOpacity = 0.2
        )
    }
  }

  m <- m |>
    addControl(
      html = paste0(
        "<b>Home range (aKDE)</b><br>Shaded regions show the estimated ",
        "home range area<br>Dashed lines show the 95% confidence interval"
      ),
      position = "topright"
    ) |>
    addLayersControl(
      baseGroups = c("World Topo Map", "World Imagery"),
      overlayGroups = ids,
      options = layersControlOptions(collapsed = FALSE)
    ) |>
    addLegend(
      position = "bottomright",
      colors = unname(cols),
      labels = ids,
      opacity = 1
    ) |>
    addEasyButton(easyButton(
      icon = "<span>Zoom to extent</span>",
      title = "Zoom to all home ranges",
      onClick = JS(
        sprintf(
          "function(btn, map) { map.fitBounds([[%f, %f], [%f, %f]]); }",
          bbox[2], bbox[1], bbox[4], bbox[3]
        )
      )
    )) |>
    # Larger text in the panels than leaflet's default (12px for the layer
    # menu, 14px for the title and legend boxes), for readability. Leaflet's
    # buttons are fixed squares, so the text button is widened to fit its label.
    htmlwidgets::prependContent(htmltools::tags$style(
      ".leaflet-control-layers, .leaflet-control.info { font-size: 16px; line-height: 1.4; }",
      ".easy-button-button { width: auto !important; padding: 0 8px; font-size: 14px; }"
    ))

  # saveWidget() can leave a folder of JavaScript libraries next to the file,
  # even though they are embedded in it. Save to a temporary folder so that
  # only the HTML file ends up with the artefacts.
  tmp_file <- file.path(tempfile("map"), "akde_map.html")
  dir.create(dirname(tmp_file))
  htmlwidgets::saveWidget(m, tmp_file, selfcontained = TRUE, title = "aKDE home ranges")
  file.copy(tmp_file, file, overwrite = TRUE)
  unlink(dirname(tmp_file), recursive = TRUE)
}
