analyzeWormz <- function(imageFolder, outputFolder, conditionsMapFile, brightfield_channel, fluorescent_channel,
                         troubleshootMode = FALSE,  # if true additional intermediate images are saved to troubleshoot
                         lightBackground = TRUE,    # findWormz() parameters...
                         threshold = "auto",
                         thresholdAdjust = 1.18,
                         fillNumPix = 0,
                         cleanNumPix = 0,
                         keepNumPix = 2000,
                         worminess_min_thr = 1.35,
                         worminess_max_thr = 2.05,
                         blurSigma = 2,
                         backgroundCorrect = TRUE,
                         showPlots = FALSE) {
  
  library(imager)
  library(tibble)
  library(imagerExtra)
  library(tools)
  library(tiff)
  library(stringr)
  library(dplyr)

  # write the set of parameters in the function call for the record
  params_df <- data.frame(lapply(as.list(match.call())[-1], eval))
  timestamp <- format(Sys.time(),"%Y%m%dT%H%M%S")
  write.csv(t(params_df), file = file.path(outputFolder, paste0("run_", timestamp, "_params.csv")))
  
  # worm names will look like cx...x_y...y_ch0z.tif
  # where x...x, y...y, and z are digits
  # x...x represents the condition
  # y...y is the image number
  # z is the channel (visible light = 0, fluorescent = 1)
  
  filenames <- Sys.glob(file.path(imageFolder, "*.tif"))
  filenames <- filenames[str_which(filenames, brightfield_channel)]  # grab just the brightfield channels
  
  #create a dataframe with analysis output
  names_analysis_output<- c("filename", "worm_number","condition_code",
            "fluo_worm",
            "fluo_background",
            "fluo_adjusted",
            "flu_sd",
            "worm_size",
            "worm_perimeter",
            "worminess")
  analysis_output <- data.frame(matrix(ncol = length(names_analysis_output), nrow = 0))
  colnames(analysis_output) <- names_analysis_output
  
  
  #loop through each ch00 file
  for (fname in filenames) {
    
    cat("\nAnalyzing", fname, "\n")
    results <- findWormz(fname,
                         lightBackground,
                         threshold,
                         thresholdAdjust,
                         fillNumPix,
                         cleanNumPix,
                         keepNumPix,
                         worminess_min_thr,
                         worminess_max_thr,
                         blurSigma,
                         backgroundCorrect,
                         showPlots)  
    
    # results list from findWormz:
    #      * final.worms is a list of masks of the pixels, one for each worm candidate
    #      * worm_imgs a list of worm images taken along the way (5 total)
    #      * stats = dataframe of statistics with one row for each final.worm in the image

    # only do stuff if worms are found in this image.  Note that findWormz will issue a warning message if none found
    if (!is.null(results) && !is.null(results$final.worms) && length(results$final.worms) > 0) {
      
      # this stuff is used to create various output filenames
      fname_basename <- basename(fname)
      fname_root     <- file_path_sans_ext(fname_basename)
      fname_extn     <- file_ext(fname_basename)
      
      # construct filename of fluorescent image data

      #locate the position in the filename root indicating the channel is brightfield and replace with the fluorescent image designation
      fname_data_image <- str_replace(fname_basename, brightfield_channel, fluorescent_channel)
      worm_data <- readTIFF(file.path(imageFolder, fname_data_image))
      condition <- str_match(fname_root, "(C\\d+)_")[,2]  # condition looks like _C1234_ in filename
      
      all_objects <- results$worm_imgs$px_all_objects  # includes all worm candidates (including touching edge), dust, etc
      background <- !all_objects
      background_mean <- mean(worm_data[background])
      
      #loop through the individual worms from this image, and add a new row to analysis_output for each worm
      for (i in 1:length(results$final.worms)) {
        
        worm <- results$final.worms[[i]]
         
        #create a new row to add to the analysis_output dataframe
        new_worm_data <- data.frame(matrix(ncol = length(names_analysis_output), nrow = 0))
        colnames(new_worm_data) <- names_analysis_output
        new_worm_data[1,] <- c(fname_data_image, 
                          i,
                          condition,
                          mean(worm_data[worm]),
                          background_mean, # background non-worm mask,
                          mean(worm_data[worm]) - background_mean,
                          sd(worm_data[worm]),
                          results$stats$numPix[i],
                          results$stats$bndryLen[i],
                          results$stats$worminess[i])
        analysis_output <- rbind(analysis_output, new_worm_data)
      }
      
      # save some of the images to output folder for documentation and troubleshooting
      # Normally only save final bw (image 4) and brightfield with color tracing (image 5))
      if (troubleshootMode) {
        start_index = 1
      } else {
        start_index = 4
      }
      for (i in start_index:5) {
        # we save the output files as jpg, not tifs, since they are much smaller
        imager::save.image(as.cimg(results$worm_imgs[[i]]), file.path(outputFolder, paste0(fname_root, "_", i, ".", "jpg") ) )
      }
    }
  }
  
  # link condition map file (effectively add 3 columns)
  conditions_map <- read.csv(conditionsMapFile)
  names(conditions_map)[1]<-"condition_code"
  analysis_output <- left_join(analysis_output, conditions_map)
  
  #save the analysis output file as a csv in your output directory
  write.csv(analysis_output, file = file.path(outputFolder, paste0("run_", timestamp, "_analyzed_data.csv")), row.names = FALSE)
  
  cat("\nDone.\n\n")
}
