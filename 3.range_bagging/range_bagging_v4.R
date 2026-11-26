
# This script takes species distribution data + covariate data to model species niches 
# and predict their distributions over multiple conditions, using a rarefaction range 
# bagging method. It then computes various range change metrics for each species.

library( here )
library( terra )
library( stringr)
library( dplyr )
library( tidyr )
library( raster )
library( nls2 ) # brute force Weibull fit for rarefaction curves
library( geometry ) # compute RB convex hulls
library( future.apply ) # parallelization inside rb.test function
library( modEvA ) # model evaluation
library( dismo )

set.seed(67402)



# GETTING DATA -----------------------------------------------------------------------------------------------------------------------------------

# import full stacks of predictors
historical <- rast(here("2.preparing_all_data","predictor_stacks","historical.tif"))
future_26_2050 <- rast(here("2.preparing_all_data","predictor_stacks","future_26_2050.tif"))
future_26_2080 <- rast(here("2.preparing_all_data","predictor_stacks","future_26_2080.tif"))
future_45_2050 <- rast(here("2.preparing_all_data","predictor_stacks","future_45_2050.tif"))
future_45_2080 <- rast(here("2.preparing_all_data","predictor_stacks","future_45_2080.tif"))
future_70_2050 <- rast(here("2.preparing_all_data","predictor_stacks","future_70_2050.tif"))
future_70_2080 <- rast(here("2.preparing_all_data","predictor_stacks","future_70_2080.tif"))

# defining named list with all predictor stacks
# the names of the SpatRaster objects are irrelevant, can be whatever user sets
# the script only assumes that the first stack is the calibration data and all the others are projection data
stacks <- mget(c(
  "historical",
  "future_26_2050", "future_26_2080",
  "future_45_2050", "future_45_2080",
  "future_70_2050", "future_70_2080"
))

# deciding on a reference stack, to obtain richness loss map
reference_stack <- "future_45_2080"

# first stack is the calibration dataset
calibration_name <- names(stacks)[1]

# remaining stacks are projection datasets
projection_names <- names(stacks)[-1]

# removing factors from stacks
stacks_nofactor <- lapply( stacks, function(x) x[[ !is.factor( stacks[[1]] ) ]] )

# importing spider occurrence data
occurrences <- vect(here("2.preparing_all_data","spiders_vector","occurrences.shp"))

# get number of records of each species
records_by_species <- as.data.frame( occurrences ) %>% 
  group_by( acceptedNa ) %>% 
  summarise( records = n())

# excluding everything below 4 records, since nls is unable to converge with less than 4 observations (which is practice means less than 3 points for the fitting of the rarefaction curve)
rb_species <- subset( records_by_species, records >= 4 )

# set model dimensionality / dimension of ranges to bag / order of models / number of covariates in each model / dimensionality of hulls to compute. d=1 not allowed
d <- 3 #TUNE

# set fixed model number caps for the rarefaction curves
# below this many unique point/variable combinations, models are enumerated exactly instead of resampled
cap_1D <- 10000 #TUNE
cap_d  <- 5000 #TUNE



# MODELLING NICHE DIMENSION AND EXPANSION -----------------------------------------------------------------------------------------------------------------------------------

# define dataframe to store loop output
niche_stats <- data.frame( 
  shortname = character(),
  real_niche_width_1D= numeric(), sampled_niche_width_1D= numeric(), portion_niche_captured_1D= numeric(), refined_1D = numeric(),
  real_niche_width_d= numeric(), sampled_niche_width_d= numeric(), portion_niche_captured_d= numeric(), inflate_d= numeric(), refined_d = numeric()
)

# define empty list for rarefaction curve data of all species
rarefaction_plot_data_1D <- list()
rarefaction_plot_data_d <- list()

# creating directories to store output files
dir.create(here("3.range_bagging","output","rarefaction_curves_1D"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("3.range_bagging","output","rarefaction_curves_d"), recursive = TRUE, showWarnings = FALSE)

# loop by species
for (w in rb_species$acceptedNa ) { 
  
  # getting shortname for folder naming
  if ( word(w, 2) == "sp." ){ # "sp" species, where numeration (third word) differentiates them (ex: Barusia sp 1, Barusia sp 2, etc)
    shortname <- paste0( word(w, 1), "_sp_", word(w, 3) ) # genus_sp_numeration
  } else {
    shortname <- paste0( word(w, 1), "_", word(w, 2) ) # genus_species
  }
  
  # subsetting all-species vector dataset to myspecies locations
  # when I exported the vector, "acceptedNameUsage" got truncated to "acceptedNa"
  myspecies_locations <- subset( occurrences, occurrences$acceptedNa == w )
  n_before <- nrow(myspecies_locations)
  
  # subsampling locations
  myspecies_locations <- spatSample( myspecies_locations, size = 1, "random", strata = stacks[[1]][[1]])
  
  n_subsampled <- nrow(myspecies_locations)
  
  # skip species which initially had sufficient locations (n_before) but not in different independent cells (n_subsampled)
  if ( n_subsampled < 4 ) { next }

  # thinning u grid to be log-spaced, to fit the rarefaction curve on: dense near the low end
  # because the rarefaction curve doesnt need a point at every single sample size, and most of the cost in a well-sampled species' curve goes to redundant points deep in the flat tail, where the function has already saturated
  u_grid_1D <- unique( round( exp( seq( log(2), log(n_subsampled), length.out = 50 ) ) ) )
  u_grid_d  <- unique( round( exp( seq( log(d+1), log(n_subsampled), length.out = 50 ) ) ) )
  
  # reporting
  print( "-----------------------------------------------------------------------------------------------------------------------------------" )
  print( paste0( "MODELLING NICHE DIMENSION AND EXPANSION OF ", shortname, " -- n: ", n_before, " -- n_subsampled: ", n_subsampled, " -- 1D curve points: ", length(u_grid_1D), " -- d curve points: ", length(u_grid_d)))


  ## GETTING COVARIATE DATA AT OBSERVATION LOCATIONS -----------------------------------------------------------------------------------------------------------------------------------
  # by extracting grid values for each point
  # removing ID and factors, since rb doesnt work with factors 
  
  # flagging index and factors as TRUE
  layers_to_rm <- c( TRUE, is.factor( stacks[[1]] ) )
  
  # extracting raster data at location points, removing ID and factors
  myspecies_env_calib <- terra::extract( stacks[[1]], myspecies_locations ) [ , -which(layers_to_rm) ] 
  
  # eliminating columns with NA
  # remove colums instead of rows because for many species entire columns are NA
  NA_cols_calib <- colSums(is.na( myspecies_env_calib )) != 0
  myspecies_env_calib <- myspecies_env_calib[ , !NA_cols_calib]
  
  # adding a very small amount of noise to covariates, so Qhull doesnt choke
  # tried adding this only to zero-variance covariates first, but its usefull to add it to all because some rows have same values and that also causes Qhull to choke
  # add noise to covariates randomized col by col, row by row
  x_reg <- myspecies_env_calib
  for (z in 1:ncol(myspecies_env_calib) ) { 
    noise <- abs( 1e-6 * myspecies_env_calib[, z][1] ) # defining noise to be orders of magnitude lower than covariate value (redundant in this dataset, in which variables are standardized)
    for ( t in 1:nrow(myspecies_env_calib) ) {  
      x_reg[t,z] <- myspecies_env_calib[t, z] + rnorm( 1, sd = noise ) 
    }
  }
  
  n <- dim(myspecies_env_calib)
  
  

  ## RAREFACTION CURVES WITH 1 DIMENSION -----------------------------------------------------------------------------------------------------------------------------------
  
  ### CONTRUCTING DF OF rarified_1D VALUES -----------------------------------------------------------------------------------------------------------------------------------
  # with niche range sampled by number of samples
  
  # set empty df
  rarified_1D <- data.frame( 
    n_locations = double(),
    mean_range = double() 
  )
  
  for ( row_i in seq_along(u_grid_1D) ){
    
    u <- u_grid_1D[row_i]
    
    # how many unique point and variable combinations
    total_combos <- n[2] * choose(n[1], u)
    
    # getting ranges
    if ( total_combos <= cap_1D ) { # few enough unique combinations to enumerate exactly instead of resampling
      
      # get unique point combinations
      point_sets <- combn( n[1], u, simplify = FALSE )
      
      # get all possible ranges
      ranges <- vector()
      k <- 1
      for (vars in 1:n[2]) {
        x0 <- x_reg[, vars]
        for (ps in point_sets) {
          x1 <- x0[ps]
          ranges[[k]] <- max(x1) - min(x1)
          k <- k + 1
        }
      }
      
    } else { # too many unique combinations, resampling preferable
      
      # get some ranges
      ranges <- vector()
      for(i in 1:cap_1D){ 
        
        # randomly subset 1 covariate 
        vars <- sample.int(n[2], size=1, replace=FALSE)
        x0 <- x_reg[,vars]
        
        # randomly subset u observations
        x1 <- x0[ sample( n[1], u, replace=FALSE ) ]
        
        # return range
        ranges[[i]] <- max(x1) - min(x1)
      }
    }
    
    # averaging ranges across all models
    mean_range <- mean(ranges)
    
    # add values to df
    rarified_1D[ row_i, 1 ] <-  u
    rarified_1D[ row_i, 2 ] <- mean_range 
  }
  
  

  ### FITTING A WEIBULL EXPONENTIAL MODEL -----------------------------------------------------------------------------------------------------------------------------------

  # define grid of potential starting values
  start_grid <- expand.grid(
    a = seq(max(rarified_1D$mean_range), # important for this specific method: the values we set to constrain "a" will also be an upper limit of the extrapolated asymptote. They should be high enough to allow good extrapolation of the asymptote. 
            max(rarified_1D$mean_range) * 4,
            length.out = 50), 
    b = seq(0.0001, 0.2, length.out = 20),
    c = seq(0.0001, 2.5, length.out = 20)) #TUNE
  
  # searching good starting values with nls2 (from discrete values on grid)
  exp_fit <- nls2( mean_range ~ a * (1 - exp(-(b * n_locations)^c)),
                   data = rarified_1D,
                   start = start_grid,
                   algorithm = "brute-force")
  
  refined_1D <- TRUE
  tryCatch({
    # refine with standard nls, starting from the coeficients of the nls2 brute-force model
    exp_fit <- nls(
      mean_range ~ a * (1 - exp(-(b * n_locations)^c)),
      data = rarified_1D,
      start = as.list(coef(exp_fit)),
      algorithm = "port", # to allow bounds
      lower = c(a = 0, b = 0, c = 0), # constrain a, b, c to positive values
      upper = c(a = Inf, b = Inf, c = Inf),
      control = nls.control(maxiter = 500) # max number of iterations
    )
  }, error = function(e) {
    refined_1D <<- FALSE
    print("Could not refine 1D rarefaction curve fit.")
    
    # define new grid of potential values
    refine_grid <- expand.grid(
      a = seq(max(rarified_1D$mean_range), # important for this specific method: the values we set to constrain "a" will also be an upper limit of the extrapolated asymptote. They should be high enough to allow good extrapolation of the asymptote.
              max(rarified_1D$mean_range) * 4,
              length.out = 100), 
      b = seq(0.0001, 0.2, length.out = 20),
      c = seq(0.0001, 2.5, length.out = 20)) #TUNE
    
    # searching refined values with nls2
    exp_fit <<- nls2( mean_range ~ a * (1 - exp(-(b * n_locations)^c)),
                      data = rarified_1D,
                      start = refine_grid,
                      algorithm = "brute-force")
  })
  
  

  ### RAREFACTION CURVE FINAL STEPS -----------------------------------------------------------------------------------------------------------------------------------
  
  # obtaining some important values
  real_niche_width_1D <- summary(exp_fit)$parameters[1,1]
  sampled_niche_width_1D <- max(rarified_1D$mean_range)
  portion_niche_captured_1D <- sampled_niche_width_1D / real_niche_width_1D
  
  # x values for plot
  xs <- seq( from = 1, to = n_subsampled, by = n_subsampled/100)
  
  # initialize image saving
  png( filename = here( "3.range_bagging","output","rarefaction_curves_1D", paste(shortname,".png")) )
  
  # plot 
  plot(rarified_1D, cex=0.2)
  preds <- data.frame ( n_locations = xs,
                        mean_range = predict( 
                          exp_fit, 
                          newdata= data.frame( n_locations = xs), 
                          type = "response")
  )
  lines(preds, col ='red')
  
  # add text to plot to inform on asymptote value
  mtext( paste0( "asymptote (value ± sd) = ", 
                 round( real_niche_width_1D, 2), 
                 "±", 
                 round( summary(exp_fit)$parameters[1,2], 2) ),
         side = 3 )
  
  # finish image saving
  dev.off()
  
  # adding rarefaction curve data to list with all species' data
  preds$species <- shortname
  rarefaction_plot_data_1D[[ length(rarefaction_plot_data_1D) +1 ]] <- preds
  
  # report 
  print("1D rarefaction curve done.")
  
  
  if( n_subsampled >= d+3 ) {
    
    ## RAREFACTION CURVES WITH TARGET DIMENSIONALITY -----------------------------------------------------------------------------------------------------------------------------------
    
    ### CONTRUCTING DF OF RARIFIED VALUES -----------------------------------------------------------------------------------------------------------------------------------
    
    plan(multisession)
    rarified_d <- future_lapply( u_grid_d, future.seed=TRUE, function(u) {
      
      # how many unique point and variable combinations
      total_combos <- choose(n[2], d) * choose(n[1], u)
      
      # fitting models
      if ( total_combos <= cap_d ) { # few enough unique combinations to enumerate exactly instead of resampling
        
        # get unique point and variable combinations
        var_sets <- combn(n[2], d, simplify = FALSE)
        point_sets <- combn(n[1], u, simplify = FALSE)
        
        # get all possible models
        all_dists <- vector()
        k <- 1
        for (vars in var_sets) {
          x0 <- x_reg[, vars]
          for (ps in point_sets) {
            x1 <- x0[ps, ]
            idx <- unique( as.vector( convhulln( x1, options = 'Pp' )))
            vertices <- x1[idx, ]
            center <- colMeans( vertices )
            distances <- as.matrix( dist( rbind( center, vertices ) ) )[-1, 1]
            all_dists[[k]] <- mean(distances)
            k <- k + 1
          }
        }
        
        # averaging mean centroid-vertice distances across all models
        mean_mean_dists <- mean( all_dists )
        
      } else {  # too many unique combinations, resampling preferable
        
        # get some models
        mean_dists_by_model <- sapply( 1:cap_d, function(i) {
          
          # randomly subset d covariates 
          vars <- sample.int(n[2], size=d, replace=FALSE)
          x0 <- x_reg[,vars]
          
          # randomly subset u observations
          x1 <- x0[ sample( n[1], u, replace=FALSE ), ]
          
          # get id's of observations which are the convex hull's vertices
          idx <- unique( as.vector( convhulln( x1, options='Pp' )))
          
          # get vertice coordinates
          vertices <- x1[idx,] 
          
          # get hull centroid
          center <- colMeans( vertices )
          
          # get distances of vertices to centroid
          distances <- as.matrix( dist( rbind( center, vertices ) ) )[-1, 1]
          
          # return
          mean(distances)
        })
        
        # averaging mean centroid-vertice distances across all models
        mean_mean_dists <- mean( mean_dists_by_model )
      }
      
      # return
      c( n_locations = u, mean_vertices_dist = mean_mean_dists )
    })
    
    rarified_d <- as.data.frame( do.call( rbind, rarified_d ))
    
    

    ### FITTING A WEIBULL EXPONENTIAL MODEL -----------------------------------------------------------------------------------------------------------------------------------

    # define grid of potential starting values
    start_grid <- expand.grid(
      a = seq(max(rarified_d$mean_vertices_dist), # important for this specific method: the values we set to constrain "a" will also be an upper limit of the extrapolated asymptote. They should be high enough to allow good extrapolation of the asymptote.
              max(rarified_d$mean_vertices_dist) * 4,
              length.out = 20), 
      b = seq(0.0001, 0.2, length.out = 5),
      c = seq(0.0001, 2.5, length.out = 5)) #TUNE
    
    # searching good starting values with nls2 (from discrete values on grid)
    exp_fit <- nls2( mean_vertices_dist ~ a * (1 - exp(-(b * n_locations)^c)),
                     data = rarified_d,
                     start = start_grid,
                     algorithm = "brute-force")
    
    refined_d <- TRUE
    tryCatch({
      # refine with standard nls, starting from the coefficients of the nls2 brute-force model
      exp_fit <- nls(
        mean_vertices_dist ~ a * (1 - exp(-(b * n_locations)^c)),
        data = rarified_d,
        start = as.list(coef(exp_fit)),
        algorithm = "port", # to allow bounds
        lower = c(a = 0, b = 0, c = 0), # constrain a, b, c to positive values
        upper = c(a = Inf, b = Inf, c = Inf),
        control = nls.control(maxiter = 500) # max number of iterations
      )
    }, error = function(e) {
      refined_d <<- FALSE
      print( paste0( "Could not refine ", d, "D rarefaction curve fit. Refining parameters using brute-force grid instead." ))
      
      # define new grid of potential values
      refine_grid <- expand.grid(
        a = seq(max(rarified_d$mean_vertices_dist), # important for this specific method: the values we set to constrain "a" will also be an upper limit of the extrapolated asymptote. They should be high enough to allow good extrapolation of the asymptote.
                max(rarified_d$mean_vertices_dist) * 4,
                length.out = 100), 
        b = seq(0.0001, 0.2, length.out = 20),
        c = seq(0.0001, 2.5, length.out = 20)) #TUNE
      
      # searching refined values with nls2
      exp_fit <<- nls2( mean_vertices_dist ~ a * (1 - exp(-(b * n_locations)^c)),
                       data = rarified_d,
                       start = refine_grid,
                       algorithm = "brute-force")
    })
    
    # get example rarefaction curve data for figure
    if (shortname == "Amaurobius_ferox") {
      save(rarified_d, file = here("3.range_bagging","output","example_rarefaction.Rdata"))
      save(exp_fit,    file = here("3.range_bagging","output","example_weibull_fit.Rdata"))
    }
    if (shortname == "Folkia_inermis") {
      save(rarified_d, file = here("3.range_bagging","output","example_rarefaction2.Rdata"))
      save(exp_fit,    file = here("3.range_bagging","output","example_weibull_fit2.Rdata"))
    }
    
    

    ### RAREFACTION CURVE FINAL STEPS -----------------------------------------------------------------------------------------------------------------------------------
    
    # obtaining some important values
    real_niche_width_d <- summary(exp_fit)$parameters[1,1]
    sampled_niche_width_d <- max(rarified_d$mean_vertices_dist)
    portion_niche_captured_d <- sampled_niche_width_d / real_niche_width_d
    
    # tolerance definition based on estimated portion of niche captured by sampling
    # the inflate value is later used to transform sampled environmental space into real (expected) environmental space
    # this means the inflate value should capture the untransformed relationship between sampled niche and real niche
    # which are approximately the max sampled niche width value and the extrapolated asymptote, respectively
    # using max() to floor inflation at 1 (cant be lower than 1)
    inflate_d <- max( real_niche_width_d / sampled_niche_width_d, 1 )
    
    # x values for plot
    xs <- seq( from = 1, to = n_subsampled, by = n_subsampled/100)
    
    # initialize image saving
    png( filename = here("3.range_bagging","output","rarefaction_curves_d", paste0(shortname, ".png")))
    
    # plot 
    plot(rarified_d, cex=0.2)
    preds <- data.frame ( n_locations = xs,
                          mean_vertices_dist = predict( exp_fit, 
                                                        newdata= data.frame( n_locations = xs), 
                                                        type = "response"))
    lines(preds, col ='red')
    
    # add text to plot to inform on asymptote value
    mtext( paste0( "asymptote (value ± sd) = ", 
                   round( real_niche_width_d, 2), 
                   "±",
                   round( summary(exp_fit)$parameters[1,2], 2) ),
           side = 3 )
    
    # finish image saving
    dev.off()
    
    # adding rarefaction curve data to list with all species' data
    preds$species <- shortname
    rarefaction_plot_data_d[[ length(rarefaction_plot_data_d) +1 ]] <- preds
    
    # report 
    print( paste0( d, "D rarefaction curve done." ))
    
    
  } else {
    real_niche_width_d <- NA
    sampled_niche_width_d <- NA
    portion_niche_captured_d <- NA
    inflate_d <- NA
    refined_d <- NA
  }
  
  
  
  ## ITERATION FINAL STEPS -----------------------------------------------------------------------------------------------------------------------------------
  
  # save results
  niche_stats[ nrow( niche_stats ) + 1, ] <- c( 
    shortname, 
    real_niche_width_1D, sampled_niche_width_1D, portion_niche_captured_1D, as.integer(refined_1D),
    real_niche_width_d, sampled_niche_width_d, portion_niche_captured_d, inflate_d, as.integer(refined_d)
    )
  
  # export output to disk
  save( niche_stats, file = here("3.range_bagging","output","niche_stats.RData") )
  save( rarefaction_plot_data_1D, file = here("3.range_bagging","output","rarefaction_plot_data_1D.RData") )
  save( rarefaction_plot_data_d, file = here("3.range_bagging","output","rarefaction_plot_data_d.RData") )
}

# reload if needed
load(here("3.range_bagging","output","niche_stats.RData"))

# convert all columns except species to numeric
niche_stats[,-1] <- lapply(niche_stats[,-1], as.numeric)

# filter out low-record species from niche_stats dataset, to fit lm
##lm_species <- subset( records_by_species, records >= 10 )
##lm_species$shortname <- sapply( lm_species$acceptedNa, function(w) {
##  if ( word(w, 2) == "sp." ){ # "sp" species, where numeration (third word) differentiates them (ex: Barusia sp 1, Barusia sp 2, etc)
##    shortname <- paste0( word(w, 1), "_sp_", word(w, 3) ) # genus_sp_numeration
##  } else {
##    shortname <- paste0( word(w, 1), "_", word(w, 2) ) # genus_species
##  }
##  shortname
##})
##niche_stats_reduced <- niche_stats[ niche_stats$shortname %in% lm_species$shortname, ]

# model
model <- lm( 
  inflate_d ~ portion_niche_captured_1D + real_niche_width_1D + sampled_niche_width_1D, 
  data = na.omit( niche_stats )
)
summary(model)

# predict target d inflate
# flooring extrapolated inflation values at 1 (cant be lower than 1)
niche_stats$inflate_d_predicted <- pmax( predict( model, niche_stats[ , c( "real_niche_width_1D", "sampled_niche_width_1D", "portion_niche_captured_1D" )]), 1 )

# check fit
with( na.omit(niche_stats), cor(inflate_d, inflate_d_predicted))



# MODELLING HABITAT SUITABILITY AND DISTRIBUTION -----------------------------------------------------------------------------------------------------------------------------------

# defining empty dataframe
# this structure is flexible to whatever raster stacks were inputted at the beginning
base_cols <- c(
  "species",
  "n",
  "n_subsampled",
  "range_size",
  "inflate",
  "boyce",
  "threshold",
  "mean_mess",
  paste0("mean_suit_", calibration_name),
  paste0("suit_pixels_", calibration_name)
)

projection_cols <- c(
  paste0("mean_suit_", projection_names),
  paste0("suit_pixels_", projection_names),
  paste0("change_", projection_names),
  paste0("p_gained_", projection_names),
  paste0("p_lost_", projection_names),
  paste0("p_site_loss_", projection_names)
)

cols <- c(base_cols, projection_cols)

stats <- setNames(
  data.frame( matrix(nrow = 0, ncol = length(cols))),
  cols
)

# creating directories to store output files
dir.create(here("3.range_bagging","output","suitability_rasters"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("3.range_bagging","output","suitability_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("3.range_bagging","output","Boyce_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("3.range_bagging","output","threshold_selection"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("3.range_bagging","output","impact_rasters"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("3.range_bagging","output","MESS_rasters"), recursive = TRUE, showWarnings = FALSE)

# loop by species
for (w in rb_species$acceptedNa ) { 
  
  # getting shortname for folder naming
  if ( word(w, 2) == "sp." ){ # "sp" species, where numeration (third word) differentiates them (ex: Barusia sp 1, Barusia sp 2, etc)
    shortname <- paste0( word(w, 1), "_sp_", word(w, 3) ) # genus_sp_numeration
  } else {
    shortname <- paste0( word(w, 1), "_", word(w, 2) ) # genus_species
  }

  # subsetting all-species vector dataset to myspecies locations
  # when I exported the vector, "acceptedNameUsage" got truncated to "acceptedNa"
  myspecies_locations <- subset( occurrences, occurrences$acceptedNa == w )
  n_before <- nrow(myspecies_locations)
  
  # subsampling locations
  myspecies_locations <- spatSample( myspecies_locations, size = 1, "random", strata = stacks[[1]][[1]])

  n_subsampled <- nrow(myspecies_locations)
  
  # skip species which initially had more than 3 locations (n_before) but dont have more than 3 in different independent cells
  if (n_subsampled <= 3 ) { next }
  
  # reporting
  print( "-----------------------------------------------------------------------------------------------------------------------------------" )
  print( paste0( "MODELLING HABITAT SUITABILITY AND DISTRIBUTION OF ", shortname, " -- n: ", n_before, " -- n_subsampled: ", n_subsampled))
  
  
  
  ## GETTING COVARIATE DATA AT OBSERVATION LOCATIONS -----------------------------------------------------------------------------------------------------------------------------------
  # by extracting grid values for each point
  # removing ID and factors, since rb doesnt work with factors 
  
  # flagging index and factors as TRUE
  layers_to_rm <- c( TRUE, is.factor( stacks[[1]] ) )

  # extracting raster data at location points, removing ID and factors
  myspecies_env_calib <- terra::extract( stacks[[1]], myspecies_locations ) [ , -which(layers_to_rm) ] 
  
  # eliminating columns with NA
  # remove colums instead of rows because for many species entire columns are NA
  NA_cols_calib <- colSums(is.na( myspecies_env_calib )) != 0
  myspecies_env_calib <- myspecies_env_calib[ , !NA_cols_calib]
  
  
  
  ## GETTING PROJECTION DATA FROM CALIBRATION AREA -----------------------------------------------------------------------------------------------------------------------------------
  
  # creating convex hull with species locations
  my_hull <- hull(myspecies_locations, type="convex") 
  
  # get range size
  range_size <- expanse( my_hull , unit='m') 
  
  # set buffer so species with wider distributions also have larger buffer
  hull_buffer <- sqrt( range_size / 10 )
  
  # adding buffer around the whole hull
  my_hull_buffered <- buffer( my_hull, width = hull_buffer) 
  
  # cropping predictor stacks to calibration area
  stacks_crop <- lapply(stacks_nofactor, crop, my_hull_buffered, mask=TRUE)
  
  # aggregating rasters if they are big
  if ( ncell( stacks_crop[[1]]) > 500*500) { #TUNE >500*500 cells means a species with a distribution bigger than 500km by 500km
    
    # deciding aggregation factor so that species with wider ranges have a bigger factor
    aggreg_factor <- sqrt( ncell( stacks_crop[[1]] ) ) *5^-3 #TUNE this value can be increased to have smaller aggregation factors and results with higher resolution
    
    # aggregating each stack
    stacks_crop <- lapply( stacks_crop, terra::aggregate, 
                           fact=aggreg_factor, fun=mean )
    
    print( paste0( "Predictor rasters were aggregated since its distribution is wider than 500x500km. Aggregation factor: ", round(aggreg_factor) ))
  }
  
  # remove columns/variables also removed in training data
  stacks_crop <- lapply( stacks_crop, function(x) x[[!NA_cols_calib]])
  
  # combine all layers from all stacks into temporary SpatRaster
  all_layers <- rast(stacks_crop)
  
  # create mask that only includes pixels with values in every layer of every stack
  na_mask <- app(all_layers, function(x) {
    if (any(is.na(x))) NA else 1
  })
  
  # apply the mask to all stacks
  # to ensure all projection layers have values in exactly same pixels, 
  # so range change statistics can be calculated between them and the raster can be uniformly populated with suitability values
  # its crucial that this is only done after removing "columns/variables also removed in training data"
  # if order reversed lots of data is lost in some cases
  stacks_crop_masked <- lapply( stacks_crop, mask, mask = na_mask)
  
  # getting covariate data from projection range
  proj_data <- lapply( stacks_crop_masked, as.data.frame, na.rm = FALSE)
  
  # remove all rows with NA values
  # its crucial that this is only done after removing "columns/variables also removed in training data"
  # if order reversed lots of data is lost in some cases
  proj_data <- lapply( proj_data, na.omit)
  
  # report
  print("Projection data was prepared.")
  
  
  
  ## rb MODELLING -----------------------------------------------------------------------------------------------------------------------------------
  
  # get niche stats of species
  specific_niche_stats <- niche_stats[ niche_stats$shortname == shortname, ]
  
  # fetch inflation factor and report
  if( n_subsampled >= d+3 ) {
    inflate <- specific_niche_stats$inflate_d
    print( paste0( "Using target dimension rarefaction curve. -- Extrapolated real niche width (asymptote): ", round( specific_niche_stats$real_niche_width_d, 2 ), " -- Sampled niche width: ", round( specific_niche_stats$sampled_niche_width_d, 2 ), " -- Portion of niche captured in sample: ", round( specific_niche_stats$portion_niche_captured_d, 2 ), " -- Inflate: ", round( specific_niche_stats$inflate_d, 2 ) ))
  } else {
    inflate <- specific_niche_stats$inflate_d_predicted
    print( paste0( "Using one dimension rarefaction curve. -- Extrapolated real niche width (asymptote): ", round( specific_niche_stats$real_niche_width_1D, 2 ), " -- Sampled niche width: ", round( specific_niche_stats$sampled_niche_width_1D, 2 ), " -- Portion of niche captured in sample: ", round( specific_niche_stats$portion_niche_captured_1D, 2 ), " -- Inflate: ", round( specific_niche_stats$inflate_d_predicted, 2 ) ))
  }
    
  v <- 200 #TUNE number of votes/models for ensemble
  p <- 0.8 #TUNE "set p = 1 to maximize the chance of sampling the extreme cases that are most informative about the niche boundary", if sink populations are unlikely, p should be close to 1 
  
  rb_models <- list()
  n <- dim(myspecies_env_calib)
  
  # adding a very small amount of noise to covariates, so Qhull doesnt choke
  # tried adding this only to zero-variance covariates first, but its usefull to add it to all because some rows have same values and that also causes Qhull to choke
  # add noise to covariates randomized col by col, row by row
  x_reg <- myspecies_env_calib
  for (z in 1:ncol(myspecies_env_calib) ) { 
    noise <- abs( 1e-6 * myspecies_env_calib[, z][1] ) # defining noise to be orders of magnitude lower than covariate value (redundant in this dataset, in which variables are standardized)
    for ( t in 1:nrow(myspecies_env_calib) ) {  
      x_reg[t,z] <- myspecies_env_calib[t, z] + rnorm( 1, sd = noise ) 
    }
  }
  
  # fitting models
  for(i in 1:v){
    
    # randomly subset d covariates 
    vars <- sample.int(n[2], size=d, replace=FALSE)
    x0 <- x_reg[,vars]
    
    # randomly subset p of observations
    x1 <- x0[ sample( n[1], ceiling( p*n[1] ), replace=FALSE ), ]
    
    # get id's of observations which are the convex hull's vertices
    idx <- unique( as.vector( convhulln( x1, options='Pp' )))
    
    # get hull vertices
    vertices <- x1[idx,]
    
    # get hull centroid
    center <- colMeans( vertices )
    
    # get inflated hull by expanding vertices using inflate value
    # ver - center: vector pointing from center to vertice
    # inflate * (ver - center): rescaled vector 
    # center + inflate * (ver - center): place vector back in space to get expanded vertice
    expanded_vertices <- t( apply( vertices, 1, function(x) center + inflate * (x - center) ))
    
    # declare model 
    rb_models[[i]] <- list( vars = vars, 
                         endpoints = expanded_vertices, # endpoints are observations which are the convex hull's vertices
                         data = x1 ) # data is bag of observations/covariates used to fit the hull/model
    
    # report
    ##print( paste0 ( shortname, " model number ", i, " of ", v, " was computed." ))
  }
  
  # report
  print("RB models were computed")
  
  # toggle
  if (FALSE) { 

    ## VISUALIZATION OF A SINGLE MODEL AS AN EXAMPLE -----------------------------------------------------------------------------------------------------------------------------------
    
    library(rgl)
    
    model <- rb_models[[1]]
    
    # plot all data points
    plot3d( model$data,
            col = "blue", size = 4 )
    
    # plot original hull
    hull <- convhulln(model$data)
    triangles3d( model$data[ t(hull), ],
                 color = "blue", alpha = 0.2 )
    
    # plot expanded vertices
    points3d( model$endpoints, 
              col = "red", size = 5 )
    
    # plot expanded hull
    expanded_hull <- convhulln(model$endpoints)
    triangles3d( model$endpoints[ t(expanded_hull), ],
                 color = "lightblue", alpha = 0.2 )
  }
  
  

  ## rb.test PROJECTION -----------------------------------------------------------------------------------------------------------------------------------

  env_suitability_scores <- lapply( proj_data, function(x_new) {

    # recovering v 
    v <- length(rb_models)
    
    # parallelize
    plan(multicore)
    
    results <- future_lapply( 1:v, function(i) {

      # subsetting covariates of x_new to only those included in this model
      test_data <- as.matrix( x_new[, rb_models[[i]]$vars] )
      
      # test if new data is inside expanded hull
      tri_pts <- tsearchn(
        rb_models[[i]]$endpoints,
        delaunayn( rb_models[[i]]$endpoints ),
        test_data
      )
      
      test_pts <- !is.na(tri_pts$p[,1])
      ##table(test_pts)
      
      # return
      test_pts
    })
    
    # report
    print("One projection was done.")
    
    out <- Reduce(`+`, results)
    return(out / v)
  })
  
  # report
  print("Suitability scores were calculated.")
  
  

  ## CREATING RASTERS OF SUITABILITY PROJECTIONS -----------------------------------------------------------------------------------------------------------------------------------
  
  # get indices of raster cells for which suitability projection exists, by getting existing rownames from projection data
  ids <- as.numeric(rownames( proj_data[[1]] ))
  
  # create empty output raster (same geometry as predictors)
  r_empty <- stacks_crop[[1]][[1]]
  values(r_empty) <- NA ; names(r_empty) <- NA ; varnames(r_empty) <- NA
  
  # create rasters and assign values vectorized
  projections <- lapply( env_suitability_scores, function(scores) {
    
    # get template
    r <- r_empty
    
    # assign suitability scores to raster geometry
    values(r)[ids] <- scores
    r
  })
  
  # making list into raster stack
  projections <- rast(projections)
  names( projections ) <- names(stacks)
  
  # export suitability projections raster
  writeRaster( projections, 
               here("3.range_bagging","output","suitability_rasters", paste0("suitability_", shortname, ".tif")),
               overwrite = TRUE)
  
  # initialize image saving
  png( filename = here("3.range_bagging","output","suitability_plots",  paste0(shortname, ".png")), width = 1000, height = 1000)
  
  # plot 
  terra::plot(projections, range=c(0, max(values(projections), na.rm = TRUE)))
  
  # finish image saving
  dev.off()
  
  # report
  print("Suitability projections plot and raster done.")
  
  

  ## MODEL EVALUATION -----------------------------------------------------------------------------------------------------------------------------------
  
  # get raster cell numbers where locations fall
  pt_cells <- cellFromXY( projections[[1]], crds(myspecies_locations) )
  
  # get backround cell indexes by doing difference between all indexes and point indexes
  all_cells <- 1:ncell(projections[[1]])
  backround_cells <- setdiff(all_cells, pt_cells)
  
  # extracting suitability predictions at location (true presence) points
  location_preds <- na.omit( values(projections[[1]])[pt_cells] )
  
  # extracting suitability prediction from non-locations (true absence) points
  backround_preds <- na.omit( values(projections[[1]])[backround_cells] )


  
  ### DETERMINISTIC DOWNSAMPLING OF BACKROUNDS PREDICTIONS -----------------------------------------------------------------------------------------------------------------------------------
  
  # sort vector before downsampling, to ensure extremes are represented exactly
  # this is crucial, as it preserves the distribution (see with histograms)
  backround_preds <- sort(backround_preds)
  
  # obtain number of backround_preds
  n_orig <- length(backround_preds)
  
  # obtain target number of backround_preds
  n_new <- n_subsampled*10
  
  if (n_orig > n_new){ # n_orig bigger than ideal, downsampling should be done
    
    # get block size
    block_size <- floor(n_orig / n_new)
    
    # truncate background predictions so they divide into blocks
    preds_trunc <- backround_preds[ 1:( block_size * n_new )]
    
    # reshape into matrix where each column is a block
    preds_matrix <- matrix( preds_trunc, nrow = block_size, ncol = n_new)
    
    # mean of each block
    preds_downsampled <- colMeans(preds_matrix)
    
  } else { # n_orig already small, downsampling not needed
    
    preds_downsampled <- backround_preds
  }

  
  # organize data to compute Boyce
  obs <- c( 
    rep(1, length(location_preds)),
    rep(0, length(preds_downsampled))
    )
  pred <- c( location_preds, preds_downsampled)
  
  # initialize image saving
  png( filename = here("3.range_bagging","output","Boyce_plots", paste0( shortname, ".png")))
  
  # compute Boyce
  boyce_result <- modEvA::Boyce(obs = obs, pred = pred,
                                plot = TRUE)

  # finish image saving 
  dev.off()
  
  # report
  print( paste0( "Ensemble was evaluated using Boyce: ", round( boyce_result$Boyce, 2 )))
  
  

  ## FIND BEST THRESHOLD TO BINARIZE -----------------------------------------------------------------------------------------------------------------------------------
  
  thresholds <- seq( from=0.01, to=0.99, by=0.01)
  
  evals <- sapply( thresholds, function(t){
    
    # binarize calib suitability projection
    binarized <- projections[[1]] > t

    # extracting presence/absence predictions at location (true presence) points
    location_preds <- na.omit( values(binarized)[pt_cells] )
    
    # extracting presence/absence prediction from non-locations (true absence) points
    backround_preds <- na.omit( values(binarized)[backround_cells] )
    
    
    
    ### DETERMINISTIC DOWNSAMPLING OF BACKROUNDS PREDICTIONS -----------------------------------------------------------------------------------------------------------------------------------
    # in the rare cases where the number of backround_preds is smaller the "n_new" target, the amount of backround preds actually gets inflated to reach "n_new"
    # this is not a problem, since T/F ratio is preserved
    
    # proportion of predicted TRUE's in backround
    p_true <- mean(backround_preds)
    
    # obtain target number of backround_preds
    n_new <- n_subsampled*10
    
    # get number of TRUE's to preserve
    n_true_new <- round(n_new * p_true) 
    
    # get downsampled backround predictions
    preds_downsampled <- c(rep(TRUE, n_true_new), rep(FALSE, n_new - n_true_new))

    # get confusion matrix values
    TP <- sum(location_preds)
    FA <- sum(!location_preds)
    FP <- sum(preds_downsampled)
    TA <- sum(!preds_downsampled)
    
    # see confusion matrix
    ##matrix(
    ##  c(TP, FA,
    ##    FP, TA),
    ##  nrow = 2,
    ##  byrow = TRUE,
    ##  dimnames = list(
    ##    actual = c("presence", "absence"),
    ##    predicted = c("presence", "absence")
    ##  )
    ##)
    
    # get evaluation metrics
    sensitivity <- TP / (TP + FA)
    specificity <- TA / (TA + FP)
    TSS <- sensitivity + specificity - 1 #TUNE eval metric can be changed
    
    # return
    TSS
  })
  
  # organize result
  thresholds_evals <- data.frame( Threshold = thresholds, TSS = evals)
  
  # get threshold that maximizes eval
  thres <- thresholds_evals$Threshold[ which.max(thresholds_evals$TSS) ]
  
  # initialize image saving
  png( filename = here("3.range_bagging","output","threshold_selection", paste0( shortname, ".png")))

  # plot evaluations by thresholds
  plot(thresholds_evals)
  abline( v = thres, col="red")
  
  # finish image saving
  dev.off()
  
  # report
  print(paste0("Optimal threshold was selected: ", thres))
  
  

  ## RANGE CHANGE STATISTICS -----------------------------------------------------------------------------------------------------------------------------------
  
  # calculate mean suitability over area, for each period/scenario
  mean_env_suitability_scores <- sapply ( env_suitability_scores, mean)
  
  # get binarized predictions
  env_suitability_bin <- lapply( env_suitability_scores, function(x){ x>thres } )
  
  # pixels kept
  kept <- setNames( vector("list", length(projection_names)), projection_names ) 
  for( i in seq_along(kept) ){
    kept[[i]] <- env_suitability_bin[[1]] & env_suitability_bin[[i+1]]
  }
  
  # pixels gained
  gained <- setNames( vector("list", length(projection_names)), projection_names ) 
  for( i in seq_along(gained) ){
    gained[[i]] <- !env_suitability_bin[[1]] & env_suitability_bin[[i+1]]
  }

  # pixels lost
  lost <- setNames( vector("list", length(projection_names)), projection_names ) 
  for( i in seq_along(lost) ){
    lost[[i]] <- env_suitability_bin[[1]] & !env_suitability_bin[[i+1]]
  }
  
  # pixels never colonized
  ##never <- setNames( vector("list", length(projection_names)), projection_names ) 
  ##for( i in seq_along(never) ){
  ##  never[[i]] <- !env_suitability_bin[[1]] & !env_suitability_bin[[i+1]]
  ##}
  
  # number of suitable, kept, gained, and lost pixels for each period/scenario
  suit_pix <- sapply( env_suitability_bin, sum ) # means number of suitable pixels allowing dispersal 
  n_kept <- sapply( kept, sum ) # means number of suitable pixels with no dispersal 
  n_gained <- sapply( gained, sum )
  n_lost <- sapply( lost, sum )
  ##n_never <- sapply( never, sum )
  
  # gained, lost, and net change percentages for each period/scenario
  p_gained <- sapply( n_gained, function(x) x /suit_pix[1] )
  p_lost <- -sapply( n_lost, function(x) x /suit_pix[1] ) # means range change percentages with no dispersal, is equivalent to "sapply( n_kept, function(x) (x - suit_pix[1] )/suit_pix[1] )"
  p_change <- p_gained + p_lost # means range change percentages allowing dispersal, is equivalent to "sapply( suit_pix[-1], function(x) (x - suit_pix[1] )/suit_pix[1] )"
  
  print("Range change statistics done.")
  
  

  ## LOCATION-WISE RESULTS -----------------------------------------------------------------------------------------------------------------------------------
  
  # extracting suitability scores at location points, filtering out occasional NAs due to aggregation
  location_suits <- na.omit( terra::extract( projections, myspecies_locations ) [,-1] )
  
  # binarizing suitability scores using same threshold as before
  location_bin <- apply( location_suits, 2, function(x){ x>thres } )
  
  # calculate fractions of predicted loss (where species was present and is predicted to be absent)
  fraction_loss <- 1 - ( colSums(location_bin) / nrow(location_bin) )
  
  # remove historical nonsense fraction_loss (this value is close to zero)
  fraction_loss <- fraction_loss[-1]
  
  

  ## CREATING RASTER OF IMPACT FOR THE SPECIES -----------------------------------------------------------------------------------------------------------------------------------

  r_impact <- r_empty
  names(r_impact) <- "impact"
  
  # getting impacted cells and transforming TRUE,FALSE to 1,0
  impacted_cells <- lost[[ which(projection_names == reference_stack) ]] *1
  
  # assign impacted/non-impacted values to raster geometry
  values(r_impact)[ids] <- impacted_cells
  
  # export raster
  writeRaster( r_impact, 
               here("3.range_bagging","output","impact_rasters", paste0("impact_", shortname, ".tif")),
               overwrite = TRUE)
  
  print("Impact raster done.")
  
  

  ## CREATING RASTER OF MESS FOR THE SPECIES -----------------------------------------------------------------------------------------------------------------------------------
  # using only the extrapolation to SSP2-4.5 2080 as proxy of extrapolation in all scenarios 

  # preparing raster stack (projection data) 
  dismo_raster <- raster::stack( stacks_crop[[5]] )
  
  # making sure the stack only contains covariates that exist in calibration data
  dismo_raster<- subset(dismo_raster, colnames(myspecies_env_calib))
  
  # compute MESS 
  mess <- dismo::mess( x=dismo_raster, 
                       v=myspecies_env_calib, 
                       full=FALSE)
  
  # export MESS raster
  writeRaster( mess, 
               here("3.range_bagging","output","MESS_rasters", paste0( "MESS_", shortname, ".tif")),
               overwrite = TRUE)
  
  

  ## GET MEAN MESS OF CALIBRATION/PROJECTION AREA -----------------------------------------------------------------------------------------------------------------------------------
  
  # get raster values
  val <- getValues(mess) 
  
  # remove inf values, which I suspect are those outside de mask (NAs in predictor stack)
  val <- val[!is.infinite(val)]
  
  # get mean of all values
  m <- mean(val, na.rm=T) 
  
  # report
  print("MESS raster and MESS mean done.")
  
  

  ## ITERATION FINAL STEPS -----------------------------------------------------------------------------------------------------------------------------------

  # adding statistics to storage dataframe 
  stats[ nrow( stats ) + 1, ] <- c( 
    shortname, n_before, n_subsampled, range_size, inflate, boyce_result$Boyce, thres, m, # base
    mean_env_suitability_scores[1], suit_pix[1], # base, for calibration dataset
    mean_env_suitability_scores[-1], suit_pix[-1], # for projection datasets
    p_change, p_gained, p_lost, fraction_loss # for projection datasets
  )
  
  # export output to disk
  save( stats, file = here("3.range_bagging","output","stats.RData"))
  
  # report
  print( paste0( shortname, " was modeled using RB."))
}



# MERGING IMPACT RASTERS TO MAKE OVERALL RASTER OF IMPACT  -----------------------------------------------------------------------------------------------------------------------------------

# fetch all raster files
impact_collection <- sprc( list.files(here("3.range_bagging","output","impact_rasters"), full.names = TRUE) )

# fetch resolution of raster with finest resolution
# will usually be 1000 by 1000
resolutions <- sapply(impact_collection, function(r) prod(res(r)))
res_target <- res(impact_collection[which.min(resolutions)])

# make template grid (biggest extent, highest resolution)
r0 <- rast(ext = ext(impact_collection), 
           resolution = res_target,
           crs =  crs(impact_collection[1]) )

# resample all rasters to same grid
impact_rs <- lapply( impact_collection, function(r) {
  resample(r, r0, method="bilinear")
})

# stack them
impact_stack <- rast(impact_rs)

# cell-wise sum
# areas with more species will have inflated impact values, but thats the point
impact_sum <- app( impact_stack, sum, na.rm = TRUE,
                   filename = here("3.range_bagging","output","impact_sum.tif"))



# MERGING MESS RASTERS TO MAKE OVERALL RASTER OF MESS -----------------------------------------------------------------------------------------------------------------------------------

# fetch all raster files
MESS_collection <- sprc( list.files(here("3.range_bagging","output","MESS_rasters"), full.names = TRUE) )

# resample all rasters to same grid
MESS_rs <- lapply( MESS_collection, function(r) {
  resample(r, r0, method="bilinear")
})

# stack them
MESS_stack <- rast(MESS_rs)

# cell-wise mean
MESS_mean <- app( MESS_stack, mean, na.rm = TRUE,
                  filename = here("3.range_bagging","output","MESS_mean.tif"))

