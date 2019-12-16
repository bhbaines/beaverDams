library(raster)
library(rgdal)
library(RColorBrewer)
library(randomForest)
library(rfUtilities)
library(glcm)
library(velox)
library(rgeos)
library(sp)
library(rgrass7)

# tell rgrass7 to use sp not stars
# use_sp()

setwd("Q:/My Drive/Coursework/MEA 592 - Biogeomorphology/Project/MEA592_ProjectData")

# # ----- Specify path to GRASS GIS installation -----
# 
# grassExecutable <- "C:/OSGeo4W64/bin/grass78.bat"
#gisBase <- "C:/OSGeo4W64/apps/grass/grass78"
# 
# # ----- Specify path to data -----
# 
# # you need to change the above to where the data is and should be on your computer
# # on Windows, it will look something like:
# # dem <-  "C:/Users/gabor/OneDrive/Plocha/R_grassgis/dem.tif"
#gisDbase <- "C:/grassdata_v2/"
# # locationPath <- "C:/grassdata/3358"
# 
# 
# # ----- Create GRASS location -----
# # load createGRASSlocation function
# # source("createGRASSlocation.R")
# 
# # pick one option (here, we are using the file we have):
# 
# # A) create a new GRASS location based on georeferenced file
# #createGRASSlocation(grassExecutable = grassExecutable,
# #                    readProjectionFrom = dem,
# #                    locationPath = locationPath)
# 
# # B) create a new GRASS location with EPSG code 4326
# # createGRASSlocation(grassExecutable = grassExecutable,
# #                     EPSG = 4326,
# #                     locationPath = locationPath)
# 
# 
# # ----- Initialisation of GRASS -----
# #initGRASS(gisBase = gisBase, 
# #          gisDbase = gisDbase,
# #          location = location,
# #          mapset = test,
# #          override = TRUE)
# 
#initGRASS(gisBase = gisBase,
#           gisDbase = gisDbase,
#          location = "nc_spm_08_grass7", 
#          mapset = "user1", 
#          SG="elevation",
#          override=T)
# 
# 
# 
# execGRASS("r.in.gdal", input=dem, output="dem")
# 
# dem <- readRAST("dem", cat=FALSE) # load DEM to R
# 
# plot(dem, main = "Digital Elevation Model", col=terrain.colors(50)) # plot DEM
# 
# execGRASS("r.topidx", input = "dem", output = "twi") # calculate topographic wetness index (TWI)
# 
# execGRASS("r.slope.aspect", elevation="dem", slope="slope", aspect="aspect") # calculate slope, aspect
# 
# execGRASS("r.info", map="aspect") # show raster info
# 
# execGRASS("r.out.gdal", input="twi", output="twi.tif", format="GTiff") # export data to GeoTIFF
# execGRASS("r.out.gdal", input="slope", output="slope.tif", format="GTiff")
# execGRASS("r.out.gdal", input="aspect", output="aspect.tif", format="GTiff")
# 
# # ----- Load environmental data to R -----
# 
# # A) Directly from GRASS GIS (wrong option)
# twi2 <- raster(readRAST("twi", cat=FALSE))
# 
# # B) Load previously saved GeoTIFF files
# twi <- raster("twi.tif")
# slope <- raster("slope.tif")
# dem <- raster("dem.tif")
# aspect <- raster("aspect.tif")
# 
# # Compare twi and twi2 - Is data source same for both layers? If not, how it may impact next modeling? 


lynch_watershed <- readOGR("watershed", "watershed", stringsAsFactors = F)
test_area <- readOGR('ML_test_data/test_area', 'lower_lynch_area', stringsAsFactors = F)

ortho_lynch_2016 <- stack('imagery/naip/lynchcreek_naip_2016.tif')
ortho_lynch_2016 <- aggregate(ortho_lynch_2016, fact=2)

red_lynch_2016 <- raster('ML_test_data/lower_lynch_test_red.tif')
green_lynch_2016 <- raster('ML_test_data/lower_lynch_test_green.tif')
blue_lynch_2016 <- raster('ML_test_data/lower_lynch_test_blue.tif')

grey_lynch_2016 <-  raster('ML_test_data/lower_lynch_test_composite.tif')

#geomorph_lynch_2016 <- raster('geomorphons/geomorph_9_1.tif')

texture_lynch_2016 <- glcm(grey_lynch_2016, 16, c(5,5))
texture_lynch_2016$glcm_correlation[is.infinite(texture_lynch_2016$glcm_correlation)] <- NA

lynch_dsm <- raster('ML_test_data/lower_lynch_test_dsm.tif')
lynch_intensity <- raster('ML_test_data/lower_lynch_test_intensity.tif')
lynch_ptdens <- raster('ML_test_data/lower_lynch_test_ptdens.tif')
lynch_slope <- raster('ML_test_data/lower_lynch_test_slope.tif')
lynch_twi <- raster('ML_test_data/lower_lynch_test_twi.tif')
lynch_ndwi <- raster('ML_test_data/lower_lynch_test_ndwi.tif')
lynch_flowacc <- raster('ML_test_data/lower_lynch_test_flwacc.tif')

test_red <- crop(red_lynch_2016, test_area)
test_green <- crop(green_lynch_2016, test_area)
test_blue <- crop(blue_lynch_2016, test_area)
test_grey <- crop(lynch_slope, test_area)
test_texture <- crop(texture_lynch_2016, test_area)
test_dsm <- crop(lynch_dsm, test_area)
test_intensity <- crop(lynch_intensity, test_area)
test_ptdens <- crop(lynch_ptdens, test_area)
test_slope <- crop(lynch_slope, test_area)
test_twi <- crop(lynch_twi, test_area)
test_ndwi <- crop (lynch_ndwi, test_area)
test_flowacc <- crop (lynch_flowacc, test_area)

#------------------------------------ RF ----------------------------------------

# Random forest classification

all_test <- stack(test_red, test_green, test_blue, test_dsm, test_intensity, test_ptdens, test_slope, test_twi, test_texture$glcm_second_moment, test_texture$glcm_contrast, test_ndwi, test_flowacc)

names(all_test) <- c("red", "green", "blue", "elev", "intensity", "ptdens", "slope", "twi", "asm", "contr", "ndwi", "flowacc")

velox_all_test <- velox(all_test)

training <- readOGR("ML_test_data/wetland_training_test", "wetland_training_test", stringsAsFactors = F)
training <- spTransform(training, CRS("+proj=lcc +lat_1=36.1666666666667 +lat_2=34.3333333333333 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m
+no_defs"))
training$class <- as.factor(training$class)
training$id <- as.numeric(training$class)

# Classes - 1 (wetland), 2 (forest), 3 (grass), 4 (bare earth)
levels(training$class) <- c('wetland', 'forest', 'grass', 'bare earth')


#Extract training data pixel values as data frame.
train_df = data.frame(matrix(vector(), nrow = 0, ncol = 12))   
for (i in 1:length(unique(training[["class"]]))){
  category <- unique(training[["class"]])[i]
  catPolys <- training[training[["class"]] == category, ]
  dataSet <- velox_all_test$extract(sp=catPolys)
  dataSet <- lapply(dataSet, function(x) {cbind(x, class = (rep(category, nrow(x))))})
  df <- do.call("rbind", dataSet)
  train_df <- rbind(train_df, df)
}

names(train_df) <- c("red", "green", "blue", "elev", "intensity", "ptdens", "slope", "twi", "asm", "contr", "ndwi", "flowacc", "class")
train_df$class <- as.factor(train_df$class)
levels(train_df$class) <- c('wetland', 'forest', 'grass', 'bare earth')

#beginCluster()
#sfQuickInit(cpus=6)
###is na.exclude the right thing to do here?
rf <- randomForest(class ~ ., data = train_df, importance=T, na.action=na.exclude)
classified_rf <- predict(all_test, rf)
rf
confusion_rf <- as.table(rf$confusion)[,1:4]
accuracy(confusion_rf)
varImpPlot(rf)

col = c('blue','green','palegreen','wheat')
plot(classified_rf, col=col, legend=F, axes = F, box = F, main = "Random forest classification")
legend(x='bottomleft', legend = levels(training$class), fill = col, bty = "n", inset=c(0.17,0), xpd=T)

writeRaster(classified_rf, filename= "ML_test_data/classified_rf.tif", format="GTiff", overwrite=T)
