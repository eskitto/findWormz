# Simple approach to image segmentation for worms when background is much lighter than worm tone
# in the visible light image, and no reflections or vignetting in the image.
#
# return value is a list with three objects:
#      * final.worms is a list of masks of the pixels, one for each worm candidate
#      * worm_imgs a list of worm images taken along the way (5 total)
#      * stats = dataframe of statistics with one row for each final.worm in the image
findWormz <- function(
  filename,
  lightBackground = TRUE,  # dark worms on light background is default
  threshold       = "auto",
  thresholdAdjust = 1,
  fillNumPix      = 0, # would use 1-5 here if needed
  cleanNumPix     = 0, # would use 1-5 here if needed
  keepNumPix      = 700,   # throw away small objects
  worminess_min_thr   = 1.35,  # threshold for perimeter to area ratio where 1 = square worm
  worminess_max_thr   = 2.05,
  blurSigma       = 2.0,  # used by isoblur
  backgroundCorrect = TRUE,
  showPlots = TRUE
) {
  
  library(imager)
  library(imagerExtra)
  library(tibble)
  library(RColorBrewer)
  library(colorspace)
  library(tools)
  library(tiff)
  
  extn <- str_to_lower(file_ext(filename))
  if (extn == "tif") {
    im <- as.cimg(readTIFF(filename))
  } else if (extn == "jpg") {
    im <- load.image(filename) # handled natively by imager library
  } else {
    cat("\nWarning: ", extn, "is an unknown file type, skipping findWormz", filename, "\n\n") 
    return(NULL)
  }
  
  if (dim(im)[4] > 1) im  <- grayscale(im) # dim 4 stores color channels, only convert to grayscale if needed
  
  im.final <- im  # will be plotted at the end (possibly with colored worms)
  final.worms <- NULL # list of worm pixsets if any worms found
  
  # small "background" correction (can create artifacts if too large lambda)
  if (backgroundCorrect) {
    im <- im %>% SPE(0.0001) 
    if (showPlots) im %>% plot()
  }
  
  # blur image then compute thresholded pixset
  px <- im %>% 
    isoblur(
      sigma = blurSigma
    ) %>% 
    threshold(
      thr = threshold, 
      adjust = thresholdAdjust
    ) %>% 
    as.pixset() 
  
  if (lightBackground) px <- !px
  
  # fill then clean to fill in some holes inside worms, eliminate some specks
  if (fillNumPix > 0 || cleanNumPix > 0) {
    px <- px %>% fill(fillNumPix) %>% clean(cleanNumPix)
  }
  
  if (showPlots) plot(px)
  # compute connected components, then discard those with too few pixels
  candidates <- split_connected(px) %>% 
    purrr::discard(~ sum(.) < keepNumPix) 
  
  # discard any components touching the 4 image edges
  # note: need to be careful about coordinates since image format in imager library is 4 dimensional
  touchesEdge <- function(im) {
    return(any(im[1,,1,1], im[dim(im)[1],,1,1], im[,1,1,1], im[,dim(im)[2],1,1]))
  }
  
  if (length(candidates) > 0) {
    if (showPlots) candidates %>% parany %>% plot()   # this plots them all
    obj.stats <- tibble(
      numPix = unlist(lapply(candidates, sum)),
      bndryLen = unlist(lapply(candidates, function(x) sum(imager::boundary(x)))),
      touchesEdge = unlist(lapply(candidates, touchesEdge)),
      worminess = bndryLen / 4 / sqrt(numPix),  # worminess is 1 for a square worm :)
      include = (worminess > worminess_min_thr & worminess < worminess_max_thr) & !touchesEdge
    )
    
    final.worms <- candidates[obj.stats$include]
    final.stats <- obj.stats %>% filter(include)
    
    print(obj.stats)
    
    if (!is.null(final.worms) && length(final.worms) > 0) {
      
      if (showPlots) final.worms %>% parany %>% plot()
      N <- length(final.worms)
      cat("Number of final worm candidates:",N,"\n")
      colors <- col2rgb(brewer.pal(min(12, max(N,3)), "Paired")) / 256.0 # color scheme - only has 12
      
      # superimpose the colored final worm candidates one by one
      for (i in 1:N) {
        df <- final.worms[[i]]
        df <- imager::where(df)
        #df <- which(final.worms[[i]])
        xbar <- round(median(df$x))
        ybar <- df %>% filter(abs(df$x - xbar) <= 1) %>% pull(y) %>% median

        im.final <- im.final %>% 
          colorise(final.worms[[i]], colors[,  1 + (i %% ncol(colors))], alpha = 0.6) %>% 
          draw_text(xbar, ybar, paste0(i, " (", round(final.stats[i, "worminess"], 2), ")"), col="red", fsize=35)
      }
      if (showPlots) plot(im.final)
      # im below includes any background processing adjustments
      # px includes blur/threshold/clean/fill
      worm_imgs <- list(im_corrected = im, px_all_objects = px, px_all_candidates = parany(candidates), px_all_final = parany(final.worms), im.final = im.final)
      
      return(list(final.worms = final.worms, worm_imgs = worm_imgs, stats = obj.stats[obj.stats$include,]))
    }
  }
  cat("\nWarning:  no worms found, returning NULL for ", filename, "\n\n")
  return(NULL)
}
