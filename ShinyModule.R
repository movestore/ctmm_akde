library(shiny)
library(ctmm)
library(sf)
library(dplyr)
library(purrr)
library(leaflet)
library(zip)
library(shinycssloaders)
library(leafem)

# Default isopleth level. Used for the slider and whenever the slider has no
# value, e.g. when the App runs in a workflow without its UI being opened.
default_level <- 0.95

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

shinyModuleUserInterface <- function(id, label) {
  ns <- NS(id) ## all IDs of UI functions need to be wrapped in ns()
  tagList(
    titlePanel("Autocorrelated Kernel Density Estimate"),
    fluidRow(
      column(4,
             sliderInput(
               ns("isopleth_levels"),
               "Isopleth level:",
               min = 0.01, max = .99, value = default_level
             ),
      ), 
      column(4, 
             checkboxGroupInput(ns("what"), "What do you want to see?", c("Lower CI" = "low", "Estimate" = "est", "Upper CI" = "high"), selected = "est")
      ),
      column(4, 
             sliderInput(ns("opacity"), "Opacity", min = 0, max = 1, value = 0.5),
             div(
               style = "display: flex; flex-wrap: wrap; gap: 6px;",
               actionButton(ns("add_all"), label = "Add all animals to map", class = "btn-sm"),
               actionButton(ns("clear"), label = "Remove all animals from map", class = "btn-sm")
             )
      ) 
    ),
    hr(),
    shinycssloaders::withSpinner(leafletOutput(ns("map"), height = "600"))
  )
}

shinyModule <- function(input, output, session, data){ ## The parameter "data" is reserved for the data object passed on from the previous app
  ns <- session$ns ## all IDs of UI functions need to be wrapped in ns()
  
  hr <- akde(data[[2]], data[[1]])
  bbx <- map(hr, ~ st_as_sf(ctmm::SpatialPolygonsDataFrame.UD(.x, level.UD = 0.97)) |> 
               st_transform(4326)) |> 
    map(st_bbox)
  bbx <- do.call(rbind, bbx)
  
  # Selected isopleth level, or the default if the UI has not been opened
  level <- reactive(input$isopleth_levels %||% default_level)

  # Home-range contours at the chosen isopleth level: a list with one sf
  # object per animal, in the same order and with the same names as hr
  hr.iso <- reactive({
    map(hr, ~ ud_contours(.x, level.UD = level()))
  })

  # Artefact: home range contours as a geopackage. Written independently of
  # the map so that it is produced even if the UI is never opened.
  observe({
    hr.iso() |>
      bind_rows() |>
      sf::st_write(appArtifactPath("homerange.gpkg"), append = FALSE)
  })

  # The same contours for the map. Each contour's name ends in "low", "est"
  # or "high" (e.g. "Cilla 95% est") so that contours can be selected interactively
  hr.leaflet <- reactive({
    map(
      hr.iso(),
      function(akde_sf) {
        akde_sf$what <- sub(".* ", "", akde_sf$name)
        st_transform(akde_sf, 4326)
      }
    )
  })
  
  output$map <- renderLeaflet({
    
    
    leaflet() |> 
      addTiles(group = "OSM") |> 
      addProviderTiles(providers$Esri.WorldImagery,
                       options = providerTileOptions(noWrap = TRUE), group = "World Imagery") |>
      fitBounds(min(bbx[, 1]), min(bbx[, 2]), max(bbx[, 3]), max(bbx[, 4])) |> 
      addHomeButton(
        ext = c(min(bbx[, 1]), min(bbx[, 2]), max(bbx[, 3]), max(bbx[, 4])), 
        group = "Full extent", position = "bottomleft")
    
  })
  
  cols <- rainbow(length(hr))
  
  observe({
    # Wait for the opacity slider before drawing. Don't wait for input$what
    # as the user may deselect all options and the map still needs to update
    req(input$opacity)
    m <- leafletProxy("map") |>
      clearShapes() |> 
      clearControls()
    for (i in 1:length(hr.leaflet())) {
      if ("est" %in% input$what) {
        m <- m |> addPolygons(data = hr.leaflet()[[i]] |> 
                                filter(what == "est"), 
                              group = names(hr.leaflet())[i],
                              label = ~name,
                              color = cols[i],
                              fillColor = cols[i],
                              fillOpacity = input$opacity)
      }
      
      if ("low" %in% input$what) {
        m <- m |> addPolygons(data = hr.leaflet()[[i]] |> 
                                filter(what == "low"),
                              group = names(hr.leaflet())[i],
                              label = ~name,
                              dashArray = "3", weight = 1,
                              opacity = 1,
                              color = "white",
                              fill = FALSE)
        
      } 
      
      if ("high" %in% input$what) {
        m <- m |> addPolygons(data = hr.leaflet()[[i]] |> 
                                filter(what == "high"),
                              group = names(hr.leaflet())[i],
                              label = ~name,
                              dashArray = "3", weight = 1,
                              color = "white",
                              opacity = 1,
                              fill = FALSE)
        
      }
    }
    m |> 
      addLayersControl(
        baseGroups = c("OSM", "World Imagery"),
        overlayGroups = names(hr.leaflet()),
        options = layersControlOptions(collapsed = TRUE)
      ) |> 
      addLegend(
        position = "bottomright",
        colors = cols,
        labels  = names(hr)) 
    
  })
  # The two buttons are left out of saved settings. Otherwise restoring them
  # would count as a click and hide or show the animals when the App opens.
  setBookmarkExclude(c("add_all", "clear"))
  observeEvent(input$add_all, {
    leafletProxy("map") |>
      showGroup(names(hr.leaflet()))
  })
  observeEvent(input$clear, {
    leafletProxy("map") |>
      hideGroup(names(hr.leaflet()))
  })
  
  
  # Artefact: summary. units = FALSE reports every area in square meters so
  # that values are comparable across animals.
  observe({
    xx <- map(hr, ~ summary(.x, level.UD = level(), units = FALSE)$CI)
    xx <- do.call(rbind, xx)
    xx |> as_tibble() |> 
      mutate(id = names(hr), unit = rownames(xx), level = level()) |> 
      select(id, level, unit, low, est, high) |> 
      write.csv(appArtifactPath("akde_summary.csv"), row.names = FALSE)
  })
  
  # Artefact: tifs
  targetDirUDs <- tempfile("uds")
  dir.create(targetDirUDs)
  
  r <- lapply(names(hr), function(x) 
    writeRaster(hr[[x]], file.path(targetDirUDs, paste0(x, ".tif")), overwrite = TRUE))
  
  zip::zip(
    zipfile = appArtifactPath("uds.zip"),
    files = list.files(targetDirUDs, full.names = TRUE, pattern = "tif$"),
    mode = "cherry-pick"
  )
  
  # List elements need names for cargo agent analyzer. Note that some apps
  # may consume this output with positional matching though.
  return(reactive({
    list(model = data[[1]], uds = hr, data = data[[2]])
  }))
}
