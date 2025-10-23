# --- Load libraries ---
library(raster)
library(sf)
library(dplyr)
library(lubridate)

# --- 1. Load SST NetCDF ---
nc_path <- "data/sst.mon.mean.nc"
sst_stack <- stack(nc_path)
sst_stack_rot <- rotate(sst_stack)  # Fix longitudes to -180:180

# Parse time info from layer names
dates <- as.Date(substr(names(sst_stack_rot), 2, 11), format = "%Y.%m.%d")
names(sst_stack_rot) <- format(dates, "%Y-%m")

# --- 2. Define stock area ID lists ---
gom_ids <- c(511, 512, 513, 514, 515)
gbk_ids <- c(522, 525, 551, 552, 561, 562)
snema_ids <- c(521, 526, 533, 534, 537, 538, 539, 611, 612, 613, 614, 615, 616,
               621, 622, 623, 624, 625, 626, 627, 628, 629, 631, 632, 633, 634,
               635, 636, 637, 638, 639)

region_lookup <- data.frame(
  Id = c(gom_ids, gbk_ids, snema_ids),
  stock = c(
    rep("GOM", length(gom_ids)),
    rep("GBK", length(gbk_ids)),
    rep("SNEMA", length(snema_ids))
  )
)

# --- 3. Load and combine shapefiles ---
shapefile_dir <- "data/Shapefiles"
shp_files <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)
default_crs <- st_crs(4326)  # WGS84 fallback

shape_list <- lapply(shp_files, function(shp) {
  sf_obj <- st_read(shp, quiet = TRUE)
  sf_obj <- st_make_valid(sf_obj)
  
  # Assign CRS if missing
  if (is.na(st_crs(sf_obj))) {
    st_crs(sf_obj) <- default_crs
  }
  
  sf_obj <- st_transform(sf_obj, st_crs(sst_stack_rot))
  
  # Detect and rename the ID column
  id_col <- grep("^id$|^Id$|id_|ID", names(sf_obj), ignore.case = TRUE, value = TRUE)
  
  if (length(id_col) == 0) {
    message("Skipping ", shp, " — no usable 'Id' column.")
    return(NULL)
  }
  
  names(sf_obj)[names(sf_obj) == id_col[1]] <- "Id"
  sf_obj <- sf_obj[, c("Id", "geometry")]
  return(sf_obj)
})

# Remove invalid or skipped shapefiles
shape_list <- Filter(Negate(is.null), shape_list)
all_shapes <- do.call(rbind, shape_list)

# --- 4. Tag shapes with stock area ---
all_shapes <- left_join(all_shapes, region_lookup, by = "Id") %>%
  filter(!is.na(stock))

# --- 5. Loop through stock regions and calculate monthly SST mean ---
results <- data.frame()

for (stock_area in unique(all_shapes$stock)) {
  message("Processing: ", stock_area)
  
  region_shape <- all_shapes %>% filter(stock == stock_area)
  region_sp <- as(region_shape, "Spatial")
  
  # Crop and mask
  sst_crop <- crop(sst_stack_rot, extent(region_sp))
  sst_mask <- mask(sst_crop, region_sp)
  
  # Calculate monthly mean SST
  mean_sst_vals <- cellStats(sst_mask, stat = "mean", na.rm = TRUE)
  
  # Build tidy result
  region_df <- data.frame(
    year = year(dates),
    month = month.name[month(dates)],
    stock = stock_area,
    var.name = "SurfaceT",
    statistic = "mean",
    value = mean_sst_vals
  )
  
  results <- bind_rows(results, region_df)
}

# --- 6. Preview and save ---
head(results)
rownames(results) <- NULL
sst <- results


write.csv(sst_anom, "data/monthly_surfaceT_by_stock.csv", row.names = FALSE)



