library(terra)
library(sf)
library(dplyr)

#load DEM raster 
dem0 <- rast("/carl-data/shared/users/sanshre/DPON_2024/Stillwater/ODM_DSM/2024_01_03_ST_M_odm_dsm_modified.tif")
dem1 <- rast("/carl-data/shared/users/sanshre/DPON_2024/Stillwater/ODM_DSM/2024_04_29_ST_M_odm_dsm.tif")
dem0_resampled <- resample(dem0, dem1)
DEM <- dem1 - dem0_resampled


shapfil <- st_read("/carl-data/shared/carl/gcp_detail/shapefile_2024_stw/shape_grid_ST_02_20_2024.shp")
shapfil <- st_transform(shapfil, crs = crs(DEM, proj = TRUE))
#load orthophoto and calculate ndvi
orthophoto <- rast("/carl-data/shared/users/sanshre/DPON_2024/Stillwater/DPON-2024-04-29/M_combined_newV/odm_orthophoto/odm_orthophoto.tif")
nir <- orthophoto[[4]]
red <- orthophoto[[1]]
ndvi <- (nir-red) / (nir+red+0.000001)

# Create empty list to store results
results <- list()

# Loop through each FID
for (fid in shapfil$FID) {
  
  # Get the polygon geometry for the current FID
  plot_geom <- shapfil[shapfil$FID == fid, ]
  
  # Crop/mask NDVI to the plot
  vi_plot <- crop(ndvi, plot_geom)
  vi_plot <- mask(vi_plot, plot_geom)
  
  # Crop/mask and resample height to the same plot
  ht_plot <- crop(DEM, plot_geom)
  ht_plot <- mask(ht_plot, plot_geom)
  ht_plot <- resample(ht_plot, vi_plot)
  
  # Calculate height threshold based on average plot plant height
  # height_threshold <- mean(values(ht_plot), na.rm = TRUE) * 0.20
  height_threshold <- quantile(values(ht_plot), probs = 0.80, na.rm = TRUE)
  
  
  # Apply subsetting to only pixels with height above threshold
  vi_sub <- values(vi_plot)[values(ht_plot) > height_threshold]
  
  # Compute plot-level statistics
  vi_stats <- list(
    FID = fid,
    mean = mean(vi_sub, na.rm = TRUE),
    sd = sd(vi_sub, na.rm = TRUE),
    skewness = e1071::skewness(vi_sub, na.rm = TRUE),
    kurtosis = e1071::kurtosis(vi_sub, na.rm = TRUE)
  )
  
  # Store results
  results[[as.character(fid)]] <- vi_stats
}
results_80 <- do.call(rbind, lapply(results, as.data.frame))
results_80 <- data.frame(FID = names(results), results_80)
write.csv(results_80, "ndvi_stats_eighty.csv", row.names = FALSE)

mean(results_80$mean, na.rm = TRUE)

#merge BYD ratings and NDVI stats 
byd_rating <- read.csv("/carl-data/home/arpoude/Desktop/gabor/FIDwithRating.csv")
full_data <- merge(results_80, byd_rating, by = "FID")

cor(full_data[, c("mean", "sd", "skewness", "kurtosis", "byd")], use = "complete.obs")
library(corrplot)
corrplot(cor(full_data[, c("mean", "sd", "skewness", "kurtosis", "byd")], use = "complete.obs"),
         method = "number", type = "upper")
plot(full_data$skewness, full_data$byd)
model <- lm(byd ~ mean + skewness + kurtosis, data = full_data)
summary(model)
