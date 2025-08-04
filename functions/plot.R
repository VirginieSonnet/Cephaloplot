#' ============================= FUNCTIONS =================================
#' Functions to plot extracted data from model output from Cephalopod. 
#' 
#' * plot_observations: plot the original observations (TODO)
#' 
#' * plot_performance: plot predictive performance and cumulative variable importance per model 
#' * plot_compare_performance: plot predictive performance and cumulative variable importance per model (TODO)
#' 
#' * plot_predictors: plot the input predictors or the predictors selected (TODO)
#' * plot_compare_predictors: plot the input predictors or the predictors selected (TODO)
#' 
#' * colmat: create a bivariate colormap
#' * plot_colmat_legend: plot the bivariate colormap as a legend 
#' * plot_prediction: plot predicted values with or without a measure of uncertainty
#' * plot_anomaly: plot difference in two prediction with or without a measure of uncertainty 
#' =============================================================================



# ----- plot_observations -----

plot_latitude_profile <- function(df,lab=NULL,latitude_curves="all",latname,valname,
                                  latitude_span=0.75,latitude_legend=TRUE,latitude_step=NULL){
  # --- Function to plot a loess latitudinal profile ---
  # Inputs: 
  #       - df: tibble, observations or predictions to use for the profile 
  #       - lab: string, label for the x axis 
  #       - latitude_curves: vector, either "all" (average+seasons),"average" or a custom vector of the northern hemisphere seasons wanted (winter = 12,1,2, spring=3,4,5, summer=6,7,8, autumn=9,10,11)
  #       - latitude_span: numeric, span for the loess function in geom_smooth (the higher the smoother), only if < 1000 points 
  #       - latitude_legend: boolean, include the legend on the bottom or not 
  # Outputs: 
  #       - plat: ggplot, smoothed latitude profile (using loess < 1000 points and gam for memory efficiency otherwise)

  # 1. Plot the smoothed curves
  # prepare plot 
  plat <- ggplot() + 
    coord_flip() + 
    # latitude labels
    scale_x_continuous(breaks = seq(-100, 100, by = 50), limits = c(-90, 90),
                       labels = function(y) ifelse(y == 0, "0°", sprintf("%d°%s", abs(y), ifelse(y < 0, "S", "N")))) + 
    labs(y=lab) +
    theme_minimal() +
    theme(panel.grid.major.y = element_line(color = "grey50", linetype = "dotted", size = 0.4),  # Light, subtle grid
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),  # Transparent background
          axis.text = element_text(size = 10, color = "gray20"), # Latitude and longitude labels
          axis.ticks = element_blank(),
          axis.title.y = element_blank()) +
    theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
  
  # plot the average latitudinal profile 
  if (any(latitude_curves %in% c("all","average"))){
    plat <- plat +  
      geom_smooth(data=df,
                  aes(x=.data[[latname]],y=.data[[valname]]),
                  color="black",method="gam",span=latitude_span,se=TRUE) 
  } 
  # plot the seasons latitudinal profile 
  if (!"average" %in% latitude_curves){
    # assign a (northern) season depending on the months
    df <- df %>%
      mutate(season=case_when(month %in% c(12,1,2)~"winter",
                              month %in% c(3,4,5)~"spring",
                              month %in% c(6,7,8)~"summer",
                              month %in% c(9,10,11)~"autumn")) %>%
      # reorder seasons
      mutate(season=factor(season,levels=c("winter","spring","summer","autumn"))) 
    
    if (!"all" %in% latitude_curves){
      # if seasons are indicated, restrict to these 
      df <- df %>%
        filter(season %in% latitude_curves) 
    }
    plat <- plat + 
      geom_smooth(data=df,
                  aes(x=.data[[latname]],y=.data[[valname]],color=season),
                  se=FALSE,method="loess",span=latitude_span) + 
      scale_color_manual(
        values = c(
          "winter" = "#1f78b4",  # deep blue
          "spring" = "#33a02c",  # green
          "summer" = "#ffcc00",  # yellow (readable on white)
          "autumn" = "#e66101"   # orange/red
        )) 
  }
  
  # Add a legend (optional)
  if (latitude_legend==TRUE){
    plat <- plat + theme(legend.position = "bottom",
                         legend.title = element_blank(),
                         legend.text = element_text(size = 10))
    # if there are more than 2 seasons, plot the legend over 2 rows 
    if ("all" %in% latitude_curves | length(latitude_curves) > 2){
      plat <- plat + guides(color = guide_legend(nrow = 2, byrow = TRUE))
    }
  }
  
  # Adjust the measurement value step (optional)
  if (!is.null(latitude_step)){
    plat <- plat +
      scale_y_continuous(breaks = seq(round(min(df[valname])), round(max(df[valname])), by = latitude_step))
  }
  return(plat)
}
  

plot_observations <- function(df_obs, colors=c("#F7FF84","#FB3C44","#2D2B63"), alpha=1,
                            correct.lightness=FALSE,lim=NULL,
                            lab="Values",title=NULL,legend_position="right",
                            land=NULL,landcolor="grey85",landfill="#E8E8E8",
                            latitude_profile=FALSE,latitude_position="left",
                            latitude_curves="all",latitude_span=0.75,
                            latitude_legend=TRUE,latitude_step=NULL){
  # --- Function to plot the prediction map with uncertainty ---
  # Inputs: 
  #       - df_obs: tibble, extracted observations data with extract_observations
  #       - colors: vector, at least two colors for the colorscale of the values 
  #       - alpha: numeric, transparency of data points 
  #       - correct.lightness: boolean, correct for lightness using the chroma package
  #       - lim: vector, low and high limits for the colorbar 
  #       - lab: string, label for the colorbar
  #       - title: string, title for the plot 
  #       - land: tibble, land mask from extract_land, if not indicated, world map is used
  #       - landcolor: string, color for the borders of the land polygon
  #       - landfill: string, color for the land polygon 
  #       - latitude_profile: boolean, if TRUE, plots a latitudinal profile 
  #       - latitude_position: string, one of left or right for the position of the latitude graph 
  #       - latitude_curves: vector, either "all" (average+seasons),"average" or a custom vector of the northern hemisphere seasons wanted (winter = 12,1,2, spring=3,4,5, summer=6,7,8, autumn=9,10,11)
  #       - latitude_span: numeric, span for the loess function in geom_smooth (the higher the smoother), only if < 1000 points 
  #       - latitude_legend: boolean, include the legend on the bottom or not 
  #       - latitude_step: numeric, step between tick labels for the mesurement value axis (x axis)
  # Outputs: 
  #       - p/pall: ggplot/patchwork plot, prediction map and latitude profile
  
  # 1. Adjust the lightness of the colors if needed 
  if(correct.lightness == TRUE){
    colors <- chroma::interp_colors(length(colors), colors=colors, interp="bezier", correct.lightness=TRUE)
  }
  
  # 2. Plot observations
  
  # 2.1. Land mask
  if (is.null(land)){
    p <- ggplot() +
      geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), 
                   fill = landfill, color = landcolor, size = 0.3)  # light gray background and borders
  } else {
    p <- ggplot() +
      geom_tile(data = land, aes(x = longitude, y = latitude), fill = landfill) 
  }
  
  # 2.2. Observations
  p <- p + 
    geom_point(data = df_obs, 
               aes(x = decimallongitude, y = decimallatitude, color = measurementvalue),
               alpha=alpha) +
    scale_color_gradientn(colours=colors,limits=lim,oob = scales::oob_squish) + 
    coord_equal() +  # Ensure proper aspect ratio but does not allow for alignment with latitude plot
  
    # 2.3. Labels 
    # longitude labels
    scale_x_continuous(breaks = seq(-200, 200, by = 50), limits = c(-180, 180),
                       labels = function(x) ifelse(x == 0, "0°", sprintf("%d°%s", abs(x), ifelse(x < 0, "W", "E")))) + 
    # latitude labels
    scale_y_continuous(breaks = seq(-100, 100, by = 50), limits = c(-90, 90),
                       labels = function(y) ifelse(y == 0, "0°", sprintf("%d°%s", abs(y), ifelse(y < 0, "S", "N")))) + 
    # axis labels
    labs(title = title,
         x = "", y = "",
         color=lab) +
    
    # 2.4. Layout
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "#D0D0D0", linetype = "dotted", size = 0.4),  # Light, subtle grid
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),  # Transparent background
          axis.text = element_text(size = 10, color = "gray20"), # Latitude and longitude labels
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position=legend_position)
  
  
  # 3. Latitude profile (optional)
  if(latitude_profile==TRUE){
    plat <- plot_latitude_profile(df=df_obs,lab=lab,latname="decimallatitude",
                                  valname="measurementvalue",
                                  latitude_curves=latitude_curves,
                                  latitude_span=latitude_span,
                                  latitude_legend=latitude_legend,
                                  latitude_step=latitude_step)
    
    # Remove the coord_equal from the map to be able to align 
    p <- p + coord_cartesian() + 
      theme(legend.margin = margin(t = -30)) # put the legend closer on the xaxis title space
    
    # 4. Patchwork of the observations plot and latitudinal profile 
    if (latitude_position=="right"){
      pall <- p + plat + plot_layout(widths=c(4,1))
    } else if (latitude_position=="left"){
      pall <- plat + p + plot_layout(widths=c(1,4))
    }
    return(pall)
  } else {
    return(p)
  }
}

# ----- plot_performance -----

plot_performance <- function(perf,plot_type="separate",metrics=NULL,
                             models=NULL,threshold=TRUE,threshold_color="green",
                             threshold_size=1,dot_color="darkblue",dot_size=4,
                             title="Model Performance"){
  # --- Function to plot the model performance metrics ---
  # Inputs: 
  #       - perf: tibble, output from extract_performance
  #       - plot_type: string, one of "separate" or "combined" to have the metrics on the same graph or next to each other
  #       - metrics: vector, metrics to plot for the separate plot (the order of the vector will be the order of the facets, and the first metric will be used to reorder the models)
  #       - threshold: boolean, plot or not vertical lines corresponding to Cephalopod thresholds
  #       - threshold_color: string, color for the vertical lines
  #       - threshold_size: numeric, linewidth for the vertical lines
  #       - dot_color: string, color for the end dots of the segments
  #       - dot_size: numeric, size for the end dots of the segments   
  #       - title: add a title to the plot 
  # Outputs: 
  #       - observations: tibble, initial measurementvalue, latitude, longitude, month, year and varfactor (if exist)
  
  # --- 0. Check that the inputs are correct options
  # plot type 
  if (any(!(plot_type %in% c("separate","combined")))){
    stop('Your plot type is not supported, choose between "separate" and "combined" (default is separate)')
  }
  # metrics
  if (any(!(metrics %in% c("R2","CUM_VIP","NSD")))){
    stop('One of your metrics does not exist: choose any or all of "R2","CUM_VIP","NSD"')
  }
  # models
  if (any(!(models %in% c("GLM","GAM","MLP","BRT","RF","SVM","ENSEMBLE")))){
    stop('One of your models does not exist: choose any or all of "GLM","GAM","MLP","BRT","RF","SVM","ENSEMBLE"')
  }
  
  
  # --- 1. Select the models to include (by default, include all)
  if (!(is.null(models))) {
    perf <- filter(perf,model %in% models) 
  }
  
  # --- 2. Reorder the models (by default by R2)
  if (!(is.null(metrics))){
    # reorder by the first metrics indicated
    perf <- mutate(perf, model=reorder(model,.data[[metrics[1]]]))
  } else {
    # by default, reorder by R2
    perf <- mutate(perf, model=reorder(model,R2))
  }
  
  # --- 3. Separate plot 
  
  if(plot_type=="separate"){
    
    # --- 3.1. Select the metrics
    # long format 
    perf <- perf %>% 
      pivot_longer(1:3,names_to="metric",values_to="value")
    
    # selection
    if (length(metrics) %in% c(1,2)){
      perf <- perf %>% 
        filter(metric %in% metrics)
    } 
    
    # vline for thresholds
    vline <- tibble(
      metric = c("R2", "CUM_VIP", "NSD"),
      vline_x = c(0.25, 50, 0.5)  # different x intercept per facet
    ) 
    
    # reorder the metrics
    if (!(is.null(metrics))){
      perf <- mutate(perf,metric=factor(metric,levels=metrics))
      vline <- mutate(vline,metric=factor(metric,levels=metrics))
    } else {
      perf <- mutate(perf,metric=factor(metric,levels=c("R2","CUM_VIP","NSD")))
      vline <- mutate(vline,metric=factor(metric,levels=c("R2","CUM_VIP","NSD")))
    }
    
    # --- 3.2. Plot 
    
    p <- perf %>%
      ggplot(aes(x = model, y = value)) +
      geom_segment(aes(xend = model, y = 0, yend = value), color = "grey80")
    
    # add vline or not 
    if(threshold==TRUE){
      p <- p + geom_hline(data = vline, aes(yintercept = vline_x),color=threshold_color,linewidth=threshold_size)
    }
    
    p <- p + 
      geom_point(size = dot_size, color = dot_color) +
      facet_wrap(.~metric,scale="free_x",strip.position = "top") + 
      coord_flip() +
      labs(
        title = title,
        x = "",
        y = ""
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.major.y = element_blank(),
            strip.text = element_text(hjust = 0, face = "bold",size=16))  # left-align and bold
  } else {
    
    # --- 4. Combined plot 
    
    p <- ggplot(perf, aes(x = model, y = R2)) +
      geom_segment(aes(xend = model, y = 0, yend = R2), color = "grey80")

    # add vline or not 
    if(threshold==TRUE){
      p <- p + geom_hline(aes(yintercept = 0.25),color=threshold_color,linewidth=threshold_size)
    }
    
    p <- p + 
      geom_point(aes(size = CUM_VIP, color = NSD)) +
      coord_flip() +
      scale_size_continuous(name = "Cum. VIP (%)") +
      scale_color_viridis_c(name = "NSD") +
      labs(
        title = title,
        x = "", y = expression(R^2)
      ) +
      theme_minimal(base_size = 13)
  }
  return(p)
}

# ----- plot_compare_performance -----

# ----- plot_predictors -----

# ----- plot_compare_predictors -----

# ----- plot_prediction -----

colmat <- function(nbreaks = 3, left=c("#F7FF84","#FB3C44","#2D2B63" ),
                   right=c("grey80","grey80"),
                   correct.lightness=c(FALSE,FALSE)) {
  # --- Function to assign color codes to raster data as a matrix ---
  # modified from: 
  # https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188#file-bivarrasterplot-r
  # Inputs: 
  #       - nbreaks: int, number of squares for the x by x square
  #       - right: vector, at least two colors for the right side of the matrix (uncertainty)
  #       - left: vector, at least two colors for the left side of the matrix (values)
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
  
  # Colors for the left side (reverse order to be intuitive from low to high values)
  my.pal.1 <- classInt::findColours(my.class, rev(left))
  # Colors for the right side (reverse order to be intuitive from low to high values)
  my.pal.2 <- classInt::findColours(my.class, rev(right))
  
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

plot_colmat_legend <- function(col.matrix,ylim=c(1,3),xlim=c(1,3),
                               ylab="Value",xlab="Uncertainty"){
  # --- Function to plot the matrix of colors ---
  # Inputs: 
  #       - col.matrix: tibble, HEXCode for the x position (xbin) and y position (ybin)
  #       - ylim: vector, low and high limits for the y axis 
  #       - xlim: vector, low and high limits for the x axis
  #       - ylab: label for the yaxis of the colorbar
  #       - xlab: label for the xaxis of the colorbar 
  # Outputs: 
  #       - p: ggplot, x by x colored square 
  # TODO: add customized numbers for x and y labels 

  p <- ggplot(col.matrix, aes(xbin, ybin, fill = HEXCode)) +
    geom_tile() +
    scale_fill_identity() +
    coord_equal(expand = FALSE) +
    # Custom x-axis tick labels (0 to xlim high in 4 steps)
    scale_x_continuous(breaks = c(0.75,5,10,15),
                       labels = round(seq(xlim[1], xlim[2], length.out=4), 2)) +
    # Custom y-axis tick labels (lowlim to highlim in 4 steps)
    scale_y_continuous(breaks = c(0.75,5,10,15),
                       labels = round(seq(ylim[1], ylim[2], length.out=4),2)) + 
    theme_void() +
    theme(aspect.ratio = 1,
          axis.title = element_text(size = 10, colour = "black",hjust = 0.5, 
                                    vjust = 1),
          axis.title.y = element_text(angle = 90, hjust = 0.5),
          axis.text=element_text(),
          axis.ticks=element_line(),
          axis.ticks.length=unit(.1, "cm")) +
    # labels
    xlab(bquote(.(xlab) ~  symbol("\256"))) + 
    ylab(bquote(.(ylab) ~  symbol("\256")))
  return(p)
}


plot_prediction <- function(output, plot_value="predictedvalue",
                            plot_uncertainty="none", uncertainty="predictednsd",
                            nbreaks=16,ycolors=c("#F7FF84","#FB3C44","#2D2B63" ),
                            xcolors=c("grey80","grey80"),correct.lightness=c(FALSE,FALSE),
                            xlim=NULL,ylim=NULL,
                            xlab="Uncertainty",ylab="Values",title=NULL,
                            legend_position="right",threshold=0.5,
                            land=NULL,landcolor="grey85",landfill="#E8E8E8",
                            stippling_bin=2.5,stippling_size=0.5,stippling_color="black",
                            anomaly=FALSE,anomaly_data,anomaly_lim_sym=TRUE,
                            contour=FALSE,contour_color="black",contour_levels=5,
                            obs=FALSE, obs_data=df_obs,obs_size=2,
                            latitude_profile=FALSE,latitude_position="left"){
    # --- Function to plot the prediction map with uncertainty ---
    # Inputs: 
    #       - output: tibble, extracted Cephalopod data with extract_prediction
    #       - plot_value: string, column pointing to the prediction to plot (predictedvalue, predictedsd, predictednsd, predictedcv,custom), default is "predictedvalue"
    #       - plot_uncertainty: string, one of "none", "bivarmap", "stippling" or "anomaly"
    #       - uncertainty: string, if plot_uncertainty is not "none", name of the 
    #                      variable to use as uncertainty
    #       - nbreaks: int, if plot_uncertainty="bivarmap", number of squares for the x by y bivariate color square
    #       - ycolors: vector, at least two colors for the colorscale of the values (if plot_uncertainty="bivarmap", for the left side of the matrix)
    #       - xcolors: vector, if plot_uncertainty="bivarmap", at least two colors for the right side of the matrix (uncertainty)
    #       - correct.lightness: boolean vector, correct for lightness using the chroma package, first element is for the left colors (values), second for the right colors when plot_uncertainty="bivarmap" (https://github.com/jiho/chroma) (default is no correction)
    #       - ylim: vector, low and high limits for the y axis colorbar (values)
    #       - xlim: vector, low and high limits for the x axis colorbar (uncertainty)
    #       - ylab: string, label for the y axis of the colorbar (values) (and if latitude_profile=TRUE for the x axis of the profile)
    #       - xlab: string, label for the x axis of the colorbar (uncertainty)
    #       - threshold: numeric, if plot_uncertainty="stippling", uncertainty value above which points are stippled
    #       - title: string, title for the plot 
    #       - land: tibble, land mask from extract_land
    #       - landcolor: string, color for the borders of the land polygon
    #       - landfill: string, color for the land polygon 
    #       - stippling_bin: numeric, bin size to round latitude and longitude and keep 1 value for stippling
    #       - stippling_size: numeric, size of the stippling dots 
    #       - stippling_color: string, color of the stippling dots 
    #       - anomaly: boolean, compute an anomaly difference between two outputs 
    #       - anomaly_data: tibble, output from extract_predictions to compute the anomaly compared to the reference output 
    #       - anomaly_lim_sim: boolean, if TRUE, adjust the ylim to be symmetrical based on max values rounded to 0.1
    #       - contour: boolean, add contour lines on the plot 
    #       - contour_color: strong, color of the contour lines
    #       - contour_levels: integer, number of bins to separate the data 
    #       - obs: boolean, plot observations data on top of the prediction map
    #       - obs_data: tibble, output from extract_observations with columns decimallatitude, decimallongitude and measurementvalue
    #       - obs_size: numeric, size of the dots for the observations
    #       - latitude_profile: boolean, if TRUE, plots a latitudinal profile 
    #       - latitude_position: string, one of left or right for the position of the latitude graph 
    # Outputs: 
    #       - p: ggplot, prediction map 
  
  
  # 1. Compute the anomaly if needed 
  if (anomaly==TRUE){
    # 1.1. Anomaly and anomaly uncertainty calculation 
    output <- output %>% {
      # rename shared columns 
      suppressMessages(rename(.,"valref"=plot_value,"uncertaintyref"=uncertainty))} %>% 
      # add the anomaly dataset 
      left_join(dplyr::select(ungroup(anomaly_data),
                              latitude,longitude,.data[[plot_value]],.data[[uncertainty]])) %>% 
      # calculate anomaly
      mutate(valanomaly=.data[[plot_value]]-valref,
             uncertaintyanomaly=sqrt(.data[[uncertainty]]^2+uncertaintyref^2)) %>% 
      # select and rename columns 
      dplyr::select(latitude,longitude,valanomaly,uncertaintyanomaly) %>% 
      rename(!!plot_value := valanomaly,
             !!uncertainty := uncertaintyanomaly)
    
    # 1.2. Update the ylim based on symmetrical max
    if(anomaly_lim_sym==TRUE & is.null(ylim)){
      maxlim <- round(max(abs(output[[plot_value]]),na.rm=TRUE)/0.1)*0.1
      ylim <- c(-maxlim,maxlim)
    }
  }
  
  # 2. Adjust the lightness of the colors if needed 
  if(correct.lightness[1] == TRUE){
    ycolors <- chroma::interp_colors(length(ycolors), colors=ycolors, interp="bezier", correct.lightness=TRUE)
  }
  if(length(correct.lightness)==2 & correct.lightness[2] == TRUE){
    xcolors <- chroma::interp_colors(length(xcolors), colors=xcolors, interp="bezier", correct.lightness=TRUE)
  }
  
  # 3. Land mask or default map world
  if (is.null(land)){
    p <- ggplot() +
      geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), 
                   fill = landfill, color = landcolor, size = 0.3,alpha=0.7)  # light gray background and borders
  } else {
    p <- ggplot() +
      geom_tile(data = land, aes(x = longitude, y = latitude), fill = landfill,alpha=0.7) 
  }
  
  # 4. Modify data based on the type of plot 
  
  # 4.1. Colormap and data bins when plot_uncertainty="bivarmap"
  if (plot_uncertainty=="bivarmap"){
    
    # 4.1.1. Colormap and legend 
    col.matrix <- colmat(nbreaks=nbreaks,right=xcolors,left=ycolors,correct.lightness=correct.lightness)
    legend <- plot_colmat_legend(col.matrix,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab)
    
    # 4.1.2. Values and uncertainty bins 
    output <- output %>% 
      mutate(ybin = cut(.data[[plot_value]],
                        breaks = seq(ylim[1], ylim[2], length.out = nbreaks + 1),  # Create 15 equal-width bins
                        labels = 1:nbreaks,  # Assign quantile numbers (1 to 15)
                        include.lowest = TRUE, 
                        right = FALSE),
             xbin = cut(.data[[uncertainty]],
                        breaks = seq(xlim[1], xlim[2], length.out = nbreaks + 1),  # Create 15 equal-width bins
                        labels = 1:nbreaks,  # Assign quantile numbers (1 to 15)
                        include.lowest = TRUE, 
                        right = FALSE)) %>% # Bins are left-closed, right-open))
      # assign all the values below the lowest ylim to be in the break 1 (same for xlim)
      # assign all the values above the highest ylim to be in the break nbreak (same for xlim)
      mutate_at(vars(ybin,xbin),as.numeric) %>% 
      mutate(ybin = case_when(.data[[plot_value]] <= ylim[1] ~ 1,
                              .data[[plot_value]] >= ylim[2] ~ nbreaks,
                              TRUE ~ ybin),
             xbin = case_when(.data[[uncertainty]] <= xlim[1] ~ 1,
                              .data[[uncertainty]] >= xlim[2] ~ nbreaks,
                              TRUE ~ xbin)) %>%
      left_join(col.matrix)
    
  # 4.2. Pattern based on threshold when plot_uncertainty="stippling
  } else if (plot_uncertainty=="stippling"){
    # define a stippling grid above the threshold and select a one point per cell
    stippling_grid <- output %>% ungroup() %>% 
      filter(.data[[uncertainty]] > threshold) %>%
      mutate(pattern="selected") %>% #just to have a label for shape
      mutate(lon_bin = round(longitude / stippling_bin) * stippling_bin,
             lat_bin = round(latitude / stippling_bin) * stippling_bin) %>%
      group_by(lon_bin, lat_bin) %>%
      slice(1) %>%  # one point per bin
      ungroup()
  }
  
  # 5. Plot the prediction 
  # 5.1. Set-up layout 
  p <- p +       
    # Labels 
    # longitude labels
    scale_x_continuous(breaks = seq(-200, 200, by = 50), limits = c(-180, 180),
                       labels = function(x) ifelse(x == 0, "0°", sprintf("%d°%s", abs(x), ifelse(x < 0, "W", "E")))) + 
    # latitude labels
    scale_y_continuous(breaks = seq(-100, 100, by = 50), limits = c(-90, 90),
                       labels = function(y) ifelse(y == 0, "0°", sprintf("%d°%s", abs(y), ifelse(y < 0, "S", "N")))) + 
    # axis labels
    labs(title = title,
         x = NULL, y = NULL,
         fill=ylab) +
    # theme
    theme_minimal() +
    theme(panel.grid.major.y = element_line(color = "grey50", linetype = "dotted", size = 0.4),  # Light, subtle grid
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),  # Transparent background
          axis.text = element_text(size = 10, color = "gray20"), # Latitude and longitude labels
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  # 5.2. Add prediction depending on the type of graph 
  # 5.2.1. Regular map for plot_uncertainty="none"
  if (plot_uncertainty=="none"){
    p <- p + 
      geom_tile(data = output, aes(x = longitude, y = latitude, fill = .data[[plot_value]])) +
      scale_fill_gradientn(colours=ycolors,limits=ylim,oob = scales::oob_squish) + 
      coord_equal() +  # Ensure proper aspect ratio but does not allow for alignment with latitude plot
      theme(legend.position=legend_position) # put the legend closer on the xaxis title space
    
    # 5.2.2. Bivariate map for plot_uncertainty="bivarmap"
  } else if (plot_uncertainty=="bivarmap"){
    p <- p + 
      # plot projection value 
      geom_tile(data = output, aes(x = longitude, y = latitude, fill = HEXCode)) +
      scale_fill_identity() +
      coord_equal() + 
      theme(legend.position="none")
    
    # add plot and legend together
    if (latitude_profile==TRUE){
      p <- p + coord_cartesian() # change coordinates to alignate with plat
    }
    if (latitude_position=="left"){
      p <- p + legend + plot_layout(widths=c(4,1))
      } else if (latitude_position=="right") {
        p <- legend + p + plot_layout(widths=c(1,4))
      }
    
    # 5.2.3. Stippling map for plot_uncertainty="stippling"
  } else if (plot_uncertainty=="stippling"){
    p <- p + 
      geom_tile(data = output, aes(x = longitude, y = latitude, fill = .data[[plot_value]])) + 
      scale_fill_gradientn(colours=ycolors,limits=ylim,oob = scales::oob_squish) + 
      coord_equal() + 
      # stippling
      geom_point(data = stippling_grid,aes(x = longitude, y = latitude,shape=pattern),size = stippling_size,alpha = 0.4,color = stippling_color) +
      scale_shape_manual(values=19,labels=xlab) + 
      labs(shape="") + 
      theme(legend.position=legend_position)
  }
  
  # 5.3. Add contour lines (optional)
  if (contour==TRUE){
    p <- p + 
      geom_contour(data=output,aes(x=longitude,y=latitude,z = .data[[plot_value]]), 
                   color = contour_color,bins=contour_levels) 
  }
  
  # 5.4. Add observations (optional)
  if (obs==TRUE & plot_value=="predictedvalue" & anomaly==FALSE){
    p <- p + 
      geom_point(data = obs_data,
                 aes(x = decimallongitude, y = decimallatitude, fill = measurementvalue),
                 shape = 21,  # filled circle with border
                 color = "black",  # border color
                 size = obs_size)        # adjust as needed
  }
  
  # 6. Add a latitudinal profile (optional)
  if(latitude_profile==TRUE){
    # 6.1. Generate latitudinal profile
    plat <- plot_latitude_profile(df=output,lab=ylab,latname="latitude",
                                  valname=plot_value,
                                  latitude_curves="average")
    
    # 6.2. Remove the coord_equal from the map to be able to align 
    p <- p + coord_cartesian()
    
    # 6.3. If the legend is on the bottom, adjust its position higher
    if (legend_position=="bottom" & plot_uncertainty != "bivarmap"){
      p <- p + theme(legend.margin = margin(t = -30))
    }
      
    # 7. Patchwork of the prediction plot, (legend) and latitudinal profile 
    if (latitude_position=="right"){
      pall <- p + plat + plot_layout(widths=c(4,1))
      } else if (latitude_position=="left"){
        pall <- plat + p + plot_layout(widths=c(1,4))
      }
    return(pall)
  } else {
    return(p)
  }
}
