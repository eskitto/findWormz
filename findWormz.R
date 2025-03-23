# Simple approach to image segmentation for worms when background is much lighter than worm tone
# in the visible light image, and no reflections or vignetting in the image.
#
# return value is a list with three objects:
#      * final.worms is a list of masks of the pixels, one for each worm candidate
#      * worm_imgs is a list of intermediate images representing each step of the image processing pipeline (5 total)
#      * stats = dataframe of statistics with one row for each final.worm in the image
#
# returns NULL if no worms are found
findWormz <- function(
  filename,
  threshold       = "auto",    # imager::threshold parameter
  thresholdAdjust = 1,         # imager::threshold parameter
  fillNumPix      = 0,         # parameter for imager::fill,  would use 1-5 here if needed
  cleanNumPix     = 0,         # parameter for imager::clean, would use 1-5 here if needed
  keepNumPix      = 700,       # throw away small objects (ie too few pixels)
  worminess_min_thr   = 1.4,  # "not thin enough" threshold 
  worminess_max_thr   = 2.7,  # "too thin" threshold
  blurSigma       = 2.0,       # used by imager::isoblur
  backgroundCorrect = TRUE,    # if TRUE call imagerExtra::SPE for contrast enhancement and to improve non-uniform lighting
  fillHoles = TRUE,            # fill all small holes inside of worms
  fillHoleMaxPct = 0.10,       # only fill holes that are less than this percentage of median worm candidate area
  showPlots = TRUE             # show intermediate plots, helpful for debugging, slows down process
) {
  
  library(imager)
  library(imagerExtra)
  library(tibble)
  library(RColorBrewer)
  library(colorspace)
  library(tools)
  library(tiff)
  library(dplyr)
  
  extn <- tolower(file_ext(filename))
  if (extn == "tif") {
    im <- as.cimg(readTIFF(filename))
  } else if (extn == "jpg") {
    im <- load.image(filename) # handled natively by imager library
  } else {
    cat("\nWarning: ", extn, "is an unknown file type, skipping findWormz", filename, "\n\n") 
    return(NULL)
  }
  
  if (dim(im)[4] > 1) im  <- grayscale(im) # dim 4 stores color channels, only convert to grayscale if needed
  
  im.final <- im      # save original grayscale image to be be plotted at the end (with colored worms overlay)
  final.worms <- NULL # will store list of worm pixsets if any worms found
  
  
  # Small "background" correction (can create artifacts if too large lambda).
  # Improves contrast and non-uniform lighting in image
  if (backgroundCorrect) {
    im <- im %>% imagerExtra::SPE(0.0001) 
    if (showPlots) im %>% plot()
  }
  
  # blur image then compute thresholded pixset
  px <- im %>% 
    imager::isoblur(
      sigma = blurSigma
    ) %>% 
    imager::threshold(
      thr = threshold, 
      adjust = thresholdAdjust
    ) %>% 
    as.pixset() 
  
  px <- !px
  
  # fill then clean to fill in some holes inside worms, eliminate some specks
  if (fillNumPix > 0 || cleanNumPix > 0) {
    px <- px %>% imager::fill(fillNumPix) %>% imager::clean(cleanNumPix)
  }
  
  if (showPlots) plot(px)
  # compute connected components, then discard those with too few pixels
  candidates <- split_connected(px) %>% 
    purrr::discard(~ sum(.) < keepNumPix) 
  
    
        
  # Used below to discard any components touching the 4 image edges
  # note: need to be careful about coordinates since image format in imager library is 4 dimensional
  touchesEdge <- function(im) {
    return(any(im[1,,1,1], im[dim(im)[1],,1,1], im[,1,1,1], im[,dim(im)[2],1,1]))
  }
  
  if (length(candidates) > 0) {
    
    if (fillHoles) {
      median.worm.area <- median(unlist(lapply(candidates, sum))) # needed to measure what is a "small" hole
      # faster to fill all holes at once vs looping over individual worms
      all.candidates <- parany(candidates)
      # calculate components of worm complement (NOTE: this excludes "background" component, which is good)
      holes <- split_connected(!all.candidates)
      # Only fill holes that have small area relative to worm area before filling.
      # This avoids filling, eg, large areas when a worms head is touching its tail,
      # forming a loop.
      small.holes <- holes %>%  purrr::discard(~ sum(.)/median.worm.area > fillHoleMaxPct) 
      if (length(small.holes) > 0) {
        if (showPlots) plot(all.candidates) # before filling holes
        # will append image of all worms to holes list before combining all in subseq line
        worms.and.holes <- small.holes
        worms.and.holes[[length(worms.and.holes)+1]] <- all.candidates 
        all.worms.filled <- parany(worms.and.holes) # this line fills the holes in all worms
        candidates <- split_connected(all.worms.filled) # split filled worms image back into indiv worms
      }
    }
    if (showPlots) candidates %>% parany %>% plot()   # this plots them all
    obj.stats <- tibble(
      numPix = unlist(lapply(candidates, sum)),
      touchesEdge = unlist(lapply(candidates, touchesEdge)),
      # Boundary pixel count is used below as a fast but *very* rough 
      # approximation for perimeter (curve length)
      bndryLen = unlist(lapply(candidates, function(x) sum(imager::boundary(x)))),
      # worminess = perimeter / 4 * sqrt(area) is a thinness or spindlyness measure
      # worminess is normalized to 1 for a square worm (aside from discretization effects)
      # In continuous setting, the smallest possible worminess is circle with value sqrt(pi)/2 (approx 0.89) by the isoperimetric inequality
      worminess = bndryLen / 4 / sqrt(numPix),  
      include = (worminess > worminess_min_thr & worminess < worminess_max_thr) & !touchesEdge
    )
    
    print(obj.stats)
    final.worms <- candidates[obj.stats$include]
    final.stats <- obj.stats %>% dplyr::filter(include)
    
    if (!is.null(final.worms) && length(final.worms) > 0) {
      
      if (showPlots) final.worms %>% parany %>% plot()
      N <- length(final.worms)
      cat("Number of final worm candidates:",N,"\n")
      colors <- col2rgb(brewer.pal(min(12, max(N,3)), "Paired")) / 256.0 # color scheme - only has 12
      
      # superimpose the colored final worm candidates one by one on the original image
      for (i in 1:N) {
        df <- imager::where(final.worms[[i]])
        xbar <- round(median(df$x))
        ybar <- df %>% dplyr::filter(abs(df$x - xbar) <= 1) %>% pull(y) %>% median

        im.final <- im.final %>% 
          colorise(final.worms[[i]], colors[,  1 + (i %% ncol(colors))], alpha = 0.6) %>% 
          draw_text(xbar, ybar, paste0(i, " (", round(final.stats[i, "worminess"], 2), ")"), col="red", fsize=35)
      }
      if (showPlots) plot(im.final)
      # im below includes any background processing adjustments
      # px_all_objects is result of the blur/threshold/clean/fill image + fill in holes of final worms.  It intentionally includes all debris.
      px_all_final <- parany(final.worms)
      worm_imgs <- list(im_corrected = im, px_all_objects = px | px_all_final, px_all_candidates = parany(candidates), px_all_final = px_all_final, im.final = im.final)
      
      return(list(final.worms = final.worms, worm_imgs = worm_imgs, stats = obj.stats[obj.stats$include,]))
    }
  }
  cat("\nWarning:  no worms found, returning NULL for ", filename, "\n\n")
  return(NULL)
}
