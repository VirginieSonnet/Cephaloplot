#' ============================= FUNCTIONS =================================
#' Functions to plot extracted data from model output from Cephalopod. 
#' 
#' * plot_observations: plot the original observations (TODO)
#' * plot_performance: plot predictive performance and cumulative variable importance per model 
#' * plot_compare_performance: plot predictive performance and cumulative variable importance per model (TODO)
#' * plot_predictors: plot the input predictors or the predictors selected 
#' * plot_compare_predictors: plot the input predictors or the predictors selected (TODO)
#' * extract_predictions: retrieve the projected observations 
#' =============================================================================


# ----- plot_observations -----

# ----- plot_performance -----

# ----- plot_compare_performance -----

# ----- plot_predictors -----

# ----- plot_compare_predictors -----

# ----- plot_predictions -----

colmat <- function(nbreaks = 3, left=c("#F7FF84","#FB3C44","#2D2B63" ),
                   right=c("grey80","grey80"),
                   correct.lightness=c(FALSE,FALSE)) {
  # --- Function to assign color codes to raster data as a matrix ---
  # modified from: 
  # https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188#file-bivarrasterplot-r
  # Inputs: 
  #       - nbreaks: int, number of squares for the x by x square
  #       - right: vector, at least two colors for the right side of the matrix
  #       - left: vector, at least two colors for the right side of the matrix
  #       - correct.lightness: boolean vector, correct for lightness using the chroma package, first element is for the left colors, second for the right colors (https://github.com/jiho/chroma)
  # Outputs: 
  #       - col.matrix: tibble with the nbreaks colors (HEXCode) for the y axis (ybin)
  #                     and the x axis (xbin) of the square
  
  
  my.data <- seq(0, 1, 1/nbreaks)
  
  # Default uses terciles (Lucchesi and Wikle [2017] doi: 10.1002/sta4.150)
  my.class <- classInt::classIntervals(my.data,
                                       n = nbreaks,
                                       style = "quantile",
  )
  
  # Colors for the left side
  if (correct.lightness[1] == TRUE){
    my.pal.1 <- classInt::findColours(my.class, chroma::interp_colors(3, colors=left, interp="bezier", correct.lightness=TRUE))
  } else {
    my.pal.1 <- classInt::findColours(my.class, left) 
  }
  
  # Colors for the right side 
  if (correct.lightness[1] == TRUE){
    my.pal.2 <- classInt::findColours(my.class, chroma::interp_colors(3, colors=right, interp="bezier", correct.lightness=TRUE))
  } else {
    my.pal.2 <- classInt::findColours(my.class, right) 
  }
  
  # Matrix
  col.matrix <- matrix(nrow = nbreaks+1, ncol = nbreaks+1, NA)
  for (i in 1:(nbreaks+1)) {
    my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
    col.matrix[(nbreaks+2) - i, ] <- classInt::findColours(my.class, my.col)
  }
  
  # (addition) Remove the first row and last column which are duplicated
  col.matrix <- col.matrix[-1,-(nbreaks+1)]
  
  # Convert to table 
  col.matrix <- col.matrix %>%
    as_tibble() %>% 
    # add the ybin as row number
    mutate(ybin = row_number()) %>% 
    # extract the xbin from the colnames as integer 
    pivot_longer(cols = -ybin, names_to = "xbin", values_to = "HEXCode") %>% 
    mutate(xbin = as.integer(str_extract(xbin,"\\d+"))) %>% 
    # add a specific bin number for each couple of ybin and xbin 
    mutate("UID" = row_number())
  
  return(col.matrix)
}

#' =============================================================================

# CHANGED ybin and xbin for ybin and xbin
plot_colmat_legend <- function(col.matrix,ylim=c(1,3),xlim=c(1,3),type="grey"){
  # --- Function to plot the matrix of colors ---
  # Inputs: 
  #       - col.matrix: tibble, HEXCode for the x position (xbin) and y position (ybin)
  #       - ylim: vector, low and high limits for the y axis 
  #       - xlim: vector, low and high limits for the x axis
  #       - type: string, either "grey" or "anomaly" for the y label  
  # Outputs: 
  #       - p: ggplot, x by x colored square 
  # TO DO: adjust the labels to not be for a 16 by 16 square + to change the x and y labs
  
  p <- ggplot(col.matrix, aes(xbin, ybin, fill = HEXCode)) +
    geom_tile() +
    scale_fill_identity() +
    coord_equal(expand = FALSE) +
    # Custom x-axis tick labels (0 to sdhighlim in 4 steps)
    scale_x_continuous(breaks = c(0.75,5,10,15),
                       labels = round(seq(xlim[1], xlim[2], length.out=4), 2)) +
    # Custom y-axis tick labels (lowlim to highlim in 4 steps)
    scale_y_continuous(breaks = c(0.75,5,10,15),
                       labels = round(seq(ylim[1], ylim[2], length.out=4))) + 
    theme_void() +
    theme(aspect.ratio = 1,
          axis.title = element_text(size = 10, colour = "black",hjust = 0.5, 
                                    vjust = 1),
          axis.title.y = element_text(angle = 90, hjust = 0.5),
          axis.text=element_text(),
          axis.ticks=element_line(),
          axis.ticks.length=unit(.1, "cm")) +
    # labels
    xlab(bquote(.("Standard deviation") ~  symbol("\256")))
  
  if(type=="grey"){
    p <- p + ylab(bquote(.("Average community median grey") ~  symbol("\256")))
  } else if(type=="anomaly"){
    p <- p + ylab(bquote(.("Anomaly: day - night") ~  symbol("\256")))
  }
  print(p)
  return(p)
}

plot_prediction <- function(output, model="ensemble",col.matrix,
                            uncertainty="predictednsd",
                            ylim=NULL,xlim=NULL,title,nbreaks,land,
                            landcolor="beige"){
    # --- Function to plot the prediction map with uncertainty ---
    # Inputs: 
    #       - extracted_prediction: tibble, extracted Cephalopod data with extract_prediction
    #                 land mask 
    #       - col.matrix: tibble, HEXCode for the x and y 
    #       - ylim: vector, low and high limits for the y axis 
    #       - xlim: vector, low and high limits for the x axis
    #       - title: string, title for the plot 
    #       - nbreaks: int, number of squares for the x by x square
    #       - landcolor: string, color for the land data points 
    # Outputs: 
    #       - p: ggplot, prediction map 
    
    # 1. Assign the grey and standard deviation values to bins 
    merged_df <- output %>% 
      mutate(ybin = cut(predictedvalue,
                           breaks = seq(ylim[1], ylim[2], length.out = nbreaks + 1),  # Create 15 equal-width bins
                           labels = 1:nbreaks,  # Assign quantile numbers (1 to 15)
                           include.lowest = TRUE, 
                           right = FALSE),
             xbin = cut(predictedsd,
                         breaks = seq(xlim[1], xlim[2], length.out = nbreaks + 1),  # Create 15 equal-width bins
                         labels = 1:nbreaks,  # Assign quantile numbers (1 to 15)
                         include.lowest = TRUE, 
                         right = FALSE)) %>% # Bins are left-closed, right-open))
      # assign all the values below the lowest ylim to be in the break 1
      # assign all the values above the highest ylim to be in the break nbreak
      mutate_at(vars(ybin,xbin),as.numeric) %>% 
      mutate(ybin = case_when(predictedvalue <= ylim[1] ~ 1,
                                 predictedvalue >= ylim[2] ~ nbreaks,
                                 TRUE ~ ybin),
             xbin = ifelse(predictedsd >= xlim[2] | predictedsd <= xlim[1],nbreaks,xbin)) %>%
      left_join(col.matrix)
    
    # 2. Plot the map 
    p <- ggplot() +
      
      # Plot projection value 
      geom_tile(data = merged_df, aes(x = longitude, y = latitude, fill = HEXCode)) +
      scale_fill_identity() + 
      
      # Overlay land mask (if needed)
      geom_tile(data = land, aes(x = x, y = y), fill = landcolor, alpha = 0.7) +
      coord_equal() +  # Ensure proper aspect ratio
      labs(title = title,
           x = "", y = "",
           fill="") +
      theme_minimal() + 
      theme(text=element_text(size=15),
            legend.position="none",
            legend.text=element_text(size=10))
    print(p)
  }
  

