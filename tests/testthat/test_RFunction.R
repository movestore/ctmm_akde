test_data <- test_data("input2_datamodel.rds")

# Write artefacts to a temporary folder instead of the working directory
artifacts_dir <- file.path(tempdir(), "artifacts")
dir.create(artifacts_dir, showWarnings = FALSE)
Sys.setenv(APP_ARTIFACTS_DIR = paste0(artifacts_dir, "/"))

test_that("isopleth levels are parsed, sorted and deduplicated", {
  expect_equal(parse_levels("0.5, 0.95"), c(0.5, 0.95))
  expect_equal(parse_levels(" 0.95 ,0.5,0.5"), c(0.5, 0.95))
  expect_equal(parse_levels("0.9"), 0.9)
  expect_equal(parse_levels(0.9), 0.9)
})

test_that("invalid isopleth levels stop the App", {
  expect_error(parse_levels(""), "No isopleth level")
  expect_error(parse_levels(NULL), "No isopleth level")
  expect_error(parse_levels("95"), "\"95\"")
  expect_error(parse_levels("0.5, abc"), "\"abc\"")
  expect_error(parse_levels("0"), "\"0\"")
  expect_error(parse_levels("1"), "\"1\"")
  expect_error(parse_levels("0.5,,0.9"), "Invalid isopleth level")
})

test_that("animal IDs are turned into safe, unique raster file names", {
  expect_equal(raster_file_names(c("Cilla", "a.b")), c("Cilla.tif", "a.b.tif"))
  expect_equal(raster_file_names("M1/2019"), "M1_2019.tif")
  expect_equal(raster_file_names(c(".hidden", "")), c("_hidden.tif", "animal.tif"))
  expect_equal(raster_file_names(c("a/b", "a:b", "A_b")), c("a_b.tif", "a_b_1.tif", "A_b_2.tif"))
  expect_equal(nchar(raster_file_names(strrep("x", 300))), 54)
})

test_that("happy path", {
  actual <- rFunction(data = test_data, isopleth_levels = "0.5, 0.95")

  # Output: the input models and data, plus one UD per animal
  expect_named(actual, c("model", "uds", "data"))
  expect_identical(actual$model, test_data[[1]])
  expect_identical(actual$data, test_data[[2]])
  expect_named(actual$uds, c("Cilla", "Mvubu"))
  expect_true(all(vapply(actual$uds, inherits, logical(1), "UD")))

  # Contours: 2 animals x 2 levels x 3 estimates
  contours <- sf::st_read(file.path(artifacts_dir, "homerange.gpkg"), quiet = TRUE)
  expect_equal(nrow(contours), 12)
  expect_setequal(contours$level, c(0.5, 0.95))
  expect_setequal(contours$estimate, c("low", "est", "high"))
  expect_true(all(contours$name == paste0(contours$id, " ", contours$level * 100, "% ", contours$estimate)))

  # Areas: one row per animal and level, larger at the higher level
  areas <- read.csv(file.path(artifacts_dir, "akde_summary.csv"))
  expect_named(areas, c("id", "level", "parameter", "low", "est", "high"))
  expect_equal(nrow(areas), 4)
  expect_true(all(areas$low <= areas$est & areas$est <= areas$high))
  cilla <- areas[areas$id == "Cilla", ]
  expect_lt(cilla$est[cilla$level == 0.5], cilla$est[cilla$level == 0.95])

  expect_gt(file.size(file.path(artifacts_dir, "akde_plots.pdf")), 0)

  # Map: a single HTML file with a hover label for each contour
  map <- readLines(file.path(artifacts_dir, "akde_map.html"), warn = FALSE)
  expect_true(any(grepl("\"Cilla, 95% range\"", map, fixed = TRUE)))
  expect_true(any(grepl("\"Mvubu, 50% range\"", map, fixed = TRUE)))
  expect_false(dir.exists(file.path(artifacts_dir, "akde_map_files")))

  # Rasters: one per animal
  expect_setequal(zip::zip_list(file.path(artifacts_dir, "uds.zip"))$filename, c("Cilla.tif", "Mvubu.tif"))
})

# In this fixture, Mvubu's movement model is not range-resident (no home range)
nonresident_data <- test_data("input2_nonresident.rds")

test_that("animals without a home range are left out of results and output", {
  # ctmm warns that it returns the model instead of a UD for Mvubu
  actual <- expect_warning(
    rFunction(data = nonresident_data, isopleth_levels = "0.5, 0.95"),
    "Fit object returned"
  )

  expect_named(actual$model, "Cilla")
  expect_named(actual$uds, "Cilla")
  expect_named(actual$data, "Cilla")

  contours <- sf::st_read(file.path(artifacts_dir, "homerange.gpkg"), quiet = TRUE)
  expect_setequal(contours$id, "Cilla")
  areas <- read.csv(file.path(artifacts_dir, "akde_summary.csv"))
  expect_setequal(areas$id, "Cilla")
  expect_setequal(zip::zip_list(file.path(artifacts_dir, "uds.zip"))$filename, "Cilla.tif")
})

test_that("the App stops if no animal has a home range", {
  only_nonresident <- list(nonresident_data[[1]]["Mvubu"], nonresident_data[[2]]["Mvubu"])
  expect_error(
    suppressWarnings(rFunction(data = only_nonresident, isopleth_levels = "0.5, 0.95")),
    "Could not estimate home range"
  )
})
