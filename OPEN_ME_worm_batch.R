#install the necessary packages--this only needs to be done the first time you run the code
install.packages("imager")
install.packages("imagerExtra")
install.packages("tibble")
install.packages("RColorBrewer")
install.packages("colorspace")
install.packages("tiff")
install.packages("tools")
install.packages("stringr")
install.packages("dplyr")

#replace with file path to directory with input, output, and R code folders
setwd("filepath")

#replace with the name at the end of your brightfield channel and fluorescent channel images
brightfield_channel <- "_ch00"
fluorescent_channel <- "_ch01"

#indicate names of three subdirectories: input has all images, output will house analyzed images, R code has findWormz and analyzeWormz
rcode_folder <- file.path(".", "R code")
image_folder <- file.path(".", "input")  # all the tifs are here
output_folder <- file.path(".", "output")   # relative to wd
if (!file.exists(output_folder)) dir.create(output_folder)

#the conditions map file should be in the input folder, and called conditions_map.csv
conditions_map_filename <- file.path(image_folder, "conditions_map.csv")

#load the findWormz and analyzeWormz functions
source(file.path(rcode_folder, "findWormz.R"))
source(file.path(rcode_folder, "analyzeWormz.R"))

#run analyzeWormz
analyzeWormz(image_folder, output_folder, conditions_map_filename, brightfield_channel, fluorescent_channel,
             troubleshootMode = FALSE, # set to TRUE for additional intermediate images (good for troubleshooting)
             threshold = "auto",
             thresholdAdjust = 1.18,
             fillNumPix = 0,
             cleanNumPix = 0,
             keepNumPix = 2000,
             worminess_min_thr = 1.4, # discard shapes that are not thin enough
             worminess_max_thr = 2.7, # discard shapes that are too thin
             blurSigma = 2,
             backgroundCorrect = TRUE,
             showPlots = FALSE)  
