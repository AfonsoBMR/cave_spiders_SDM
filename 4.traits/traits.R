
# This script models the relationship between functional traits and climate change vulnerability.

library( here )
library( readxl )
library( tibble )
library( ggplot2 )
library( missForest ) # imputation
library( GGally ) 
library( stringr )
library( glmmTMB ) # GLMM
library( ggeffects ) # plot model responses
library( DHARMa ) # model residual plots
library( performance ) # model diagnostics
library( terra )
library( BAT ) # Gower distance, functional distance trees
library( phyloregion ) # PD/FD calculation
library( Matrix ) # sparse matrices

set.seed(84782)



# PROCESS VULNERABILITY DATA -----------------------------------------------------------------------------------------------------------------------------------

# import climate change vulnerability statistics
load(here::here("3.range_bagging","output","stats.RData"))

# convert all columns except species to numeric
stats[,-1] <- lapply(stats[,-1], as.numeric)

# filter by Boyce, eliminating rows whose stats are based on poor-performing models
keep <- stats$boyce > 0.5 #TUNE boyce threshold 
stats_sub <- stats[ keep, ] 
nrow(stats) ; nrow(stats_sub) 

# drop all "historical", "suitability" and "gained" statistics
stats_sub <- stats_sub[ , !( grepl( "historical|suit|gained", colnames(stats_sub))) ]

# create range "retained" variables
stats_sub$retained_26_2080 <- stats_sub$change_future_26_2080 + 1
stats_sub$retained_45_2080 <- stats_sub$change_future_45_2080 + 1
stats_sub$retained_70_2080 <- stats_sub$change_future_70_2080 + 1
# 0 -> total range loss 
# 1 -> range size stability
# >1 -> range increase
# this way the distribution is positive and zero-inflated, instead of inflated at -1

# inspect "retained" distributions
hist( stats_sub$retained_26_2080, breaks = 50 )
hist( stats_sub$retained_45_2080, breaks = 50 )
hist( stats_sub$retained_70_2080, breaks = 50 )



# PROCESS TRAITS DATA -----------------------------------------------------------------------------------------------------------------------------------

# import functional traits database
traits_all <- as.data.frame(read_excel(here::here("4.traits","Database_Mammola_et_al_2022_Figshare_V3.xlsx")))

# uniformizing species name structure
traits_all$Genus_species <- gsub( " ", "_", traits_all$Genus_species)
traits_all$Genus_species <- gsub( "\\.", "", traits_all$Genus_species)
traits_all[ traits_all$Genus_species == "Gen_sp_1", "Genus_species" ] <- "Dysderidae_Gen."
traits_all[ traits_all$Genus_species == "Gen_12_sp_1", "Genus_species" ] <- "Linyphiidae_Gen."

# assign data types to columns
traits_all[ c(13:24, 59, 66:75) ] <- lapply( traits_all[ c(13:24, 59, 66:75) ], function(x) as.factor(as.character(x)))
traits_all[ c(25:58, 60:65) ] <- lapply( traits_all[ c(25:58, 60:65) ], function(x) as.numeric(as.character(x)))

# check amount of NA's
na_count <- colSums(is.na(traits_all)) ; na_count

# drop columns with excess proportion of NA's
traits_all <- traits_all[ na_count/nrow(traits_all) < 0.33 ] #TUNE

# complete NA's with RF imputation
imp <- missForest( traits_all[, 13:55 ], maxiter = 1000 )
traits_all[, 13:55 ] = imp$ximp
    
# check imputation quality, simple method
# informs roughly how noisy imputation was
imp$OOBerror 

# check imputation quality, complex method
data_sub <- traits_all[, 13:55]
truth <- data_sub[complete.cases(data_sub), ]

num_traits <- names(truth)[sapply(truth, is.numeric)]
cat_traits <- names(truth)[sapply(truth, is.factor)]

results_num <- list()
results_cat <- list()

for (i in 1:50) {
  masked <- prodNA(truth, noNA = 0.2)
  imp_i  <- missForest(masked, xtrue = truth, verbose = FALSE)
  
  # numeric traits: Q, zRMSE, intercept/slope
  for (tr in num_traits) {
    idx <- is.na(masked[[tr]])
    if (sum(idx) < 5) next
    obs <- truth[idx, tr]
    est <- imp_i$ximp[idx, tr]
    
    fit <- lm(est ~ obs)
    Q <- 1 - sum((est - obs)^2) / sum((obs - mean(obs))^2)
    zrmse <- sqrt(mean(((est - mean(obs))/sd(obs) - (obs - mean(obs))/sd(obs))^2))
    
    results_num[[length(results_num)+1]] <- data.frame(
      iter = i, trait = tr,
      intercept = coef(fit)[1], slope = coef(fit)[2],
      Q = Q, zRMSE = zrmse
    )
  }
  
  # categorical traits: proportion correctly classified
  for (tr in cat_traits) {
    idx <- is.na(masked[[tr]])
    if (sum(idx) < 5) next
    obs <- truth[idx, tr]
    est <- imp_i$ximp[idx, tr]
    
    pcc <- mean(as.character(est) == as.character(obs))
    
    results_cat[[length(results_cat)+1]] <- data.frame(
      iter = i, trait = tr, PCC = pcc
    )
  }
}

results_num_df <- do.call(rbind, results_num)
results_cat_df <- do.call(rbind, results_cat)

# check results for continuous traits
aggregate(cbind(intercept, slope, Q, zRMSE) ~ trait, results_num_df, mean)
# intercept = 0, slope = 1 -> accurate
# Q close to 1 -> precise (Q < 0 means the trait mean is a better guess than imputation)
# zRMSE < 1 -> imputed values beat just using the mean; zRMSE > 1 means they're worse

# every trait imputs well except Sexual_Size_Dimorphism; AME, Prosoma_shape mediocre imputation, still ok for rough patterns

# results for categorical traits
pcc_summary <- aggregate(PCC ~ trait, results_cat_df, mean)

# get baseline: accuracy of just always guessing the most common class
# compare each PCC to its traits baseline: if PCC is only marginally above baseline, imputation isnt adding information
baseline <- sapply(truth[cat_traits], function(x) max(prop.table(table(x))))
baseline_df <- data.frame(trait = names(baseline), baseline = as.numeric(baseline))

# compare
print(merge(baseline_df,pcc_summary, by = "trait"))

# Ambush_hunter, Orb_web, Alien are not meaningfully imputed due to extreme class imbalance. still ok for rough patterns



# SPECIFIC TRAITS AS PREDICTORS OF VULNERABILITY -----------------------------------------------------------------------------------------------------------------------------------

## FURTHER PROCESS TRAITS DATA -----------------------------------------------------------------------------------------------------------------------------------

# keep only traits related to subterranean specialization 
traits <- traits_all[, c( 
  "Family", "Genus", "Genus_species", 
  "Ecological_classification", "Pigment", "Eyeless", "Eyes_regression", "Leg_elongation", 
  "Body_length_avg", "Prosoma_length_avg", "Prosoma_width_avg", "Femur_I_length_avg", "Femur_elongation", "Tibia_I_length_avg", "Tibia_elongation"
)]

# plot continuous traits to check distributions and colinearity
GGally::ggpairs( as.data.frame( traits[ , sapply( traits, is.numeric )]))

# checking for outliers with dotplots
par(mfrow=c(3,3)) #TUNE
for (col_name in names( traits[ , sapply( traits, is.numeric )] )) {
  dotchart(
    traits[[col_name]], 
    main = col_name
  )
}

# plotting factor traits and checking for unbalance between classes
facs <- traits[, sapply(traits, is.factor)]
par(mfrow = c(2, 3), mar = c(2, 2, 2, 2)) #TUNE
for (i in seq_along(facs)) {
  plot(facs[[i]],  main = names(facs)[i])
}

# turn "Pigment" factor into continuous score
traits$pigment_adp <- as.numeric( stringr::str_replace_all( 
  traits$Pigment,
  c("depigmented" = "1",
    "Depigmented" = "1",
    "Partly" = "0",
    "Variable" = "0",
    "Fully" = "-1",
    "Pigmented" = "-1" )
))
# -1 -> fully pigmented
# 0 -> partly pigmented or variable
# 1 -> depigmented

# concatenate "Eyeless" and "Eyes_regression" into single continuous score
traits$eye_adp <- ifelse(
  traits$Eyes_regression == 1, 0, ifelse(
    traits$Eyeless == 1, 1, -1 ))
# -1 -> normal eyes
# 0 -> regressed eyes
# 1 -> eyeless

# keep only final chosen traits
traits <- traits[, c( "Family", "Genus", "Genus_species", "Body_length_avg", "Femur_elongation", "pigment_adp", "eye_adp" )]

# col id's of numeric traits
numeric <- sapply( traits, is.numeric )

# get SD's of numerical parts of dataframe before scaling
##sds <- apply( traits[ , numeric], 2, sd, na.rm=TRUE)

# standardize numerical part of dataframe
##traits[ , numeric] <- scale( 
##  traits[ , numeric ], 
##  center = TRUE, scale = TRUE
##  )

# get SD's of "Body_length_avg", "Femur_elongation" before scaling
sds <- apply( traits[ , c("Body_length_avg", "Femur_elongation")], 2, sd, na.rm=TRUE)

# standardize only "Body_length_avg" and "Femur_elongation"
# because "pigment_adp" and "eye_adp" are already very close to being mean=0 and with SD=1, and this way their coefficients remain interpretable
traits[ , c("Body_length_avg", "Femur_elongation")] <- scale( 
  traits[ , c("Body_length_avg", "Femur_elongation")], 
  center = TRUE, scale = TRUE
  )



## MERGE DATASETS -----------------------------------------------------------------------------------------------------------------------------------

# merging both datasets
merged <- merge( 
  stats_sub, traits, 
  by.x = "species", by.y = "Genus_species",
  all.x = TRUE)

# plot retained against traits
par(mfrow=c(2,2))
with(merged, plot(retained_26_2080 ~ Body_length_avg, col="darkgreen", pch=16))
with(merged, points(retained_70_2080 ~ Body_length_avg, col="darkred", pch=16))
with(merged, plot(retained_26_2080 ~ Femur_elongation, col="darkgreen", pch=16))
with(merged, points(retained_70_2080 ~ Femur_elongation, col="darkred", pch=16))
with(merged, plot(retained_26_2080 ~ pigment_adp, col="darkgreen", pch=16))
with(merged, points(retained_70_2080 ~ pigment_adp, col="darkred", pch=16))
with(merged, plot(retained_26_2080 ~ eye_adp, col="darkgreen", pch=16))
with(merged, points(retained_70_2080 ~ eye_adp, col="darkred", pch=16))



## MODELLING -----------------------------------------------------------------------------------------------------------------------------------

# fit RCP2.6 Tweedie GLMM
m_tweedie_26 <- glmmTMB::glmmTMB(
  retained_26_2080 ~ 
    Body_length_avg + Femur_elongation + pigment_adp + eye_adp +
    (1 | Family/Genus),
  family = tweedie(link = "log"),
  data = merged
)

# fit RCP4.5 Tweedie GLMM
m_tweedie_45 <- glmmTMB::glmmTMB(
  retained_45_2080 ~ 
    Body_length_avg + Femur_elongation + pigment_adp + eye_adp +
    (1 | Family/Genus),
  family = tweedie(link = "log"),
  data = merged
)

# fit RCP7.0 Tweedie GLMM
m_tweedie_70 <- glmmTMB::glmmTMB(
  retained_70_2080 ~ 
    Body_length_avg + Femur_elongation + pigment_adp + eye_adp +
    (1 | Family/Genus),
  family = tweedie(link = "log"),
  data = merged
)

# fit respective null models (just phylogeny random effects, no trait fixed effects)
m_null_26 <- glmmTMB::glmmTMB(
  retained_26_2080 ~ 
    1 + (1|Family/Genus),
  family = tweedie(link="log"),
  data = merged
)
m_null_45 <- glmmTMB::glmmTMB(
  retained_45_2080 ~ 
    1 + (1|Family/Genus),
  family = tweedie(link="log"),
  data = merged
)
m_null_70 <- glmmTMB::glmmTMB(
  retained_70_2080 ~ 
    1 + (1|Family/Genus),
  family = tweedie(link="log"),
  data = merged
)

# compare null and full tweedie models, full model AIC should be considerably lower
AIC( m_null_26, m_tweedie_26, m_null_45, m_tweedie_45, m_null_70, m_tweedie_70 )

# inspect fits and residuals using DHARMa 
plot( DHARMa::simulateResiduals( m_tweedie_26  ))
plot( DHARMa::simulateResiduals( m_tweedie_45  ))
plot( DHARMa::simulateResiduals( m_tweedie_70 ))

# check collinearity
performance::check_collinearity( m_tweedie_26 )
performance::check_collinearity( m_tweedie_45 )
performance::check_collinearity( m_tweedie_70 )

# see original covariate SD's before scaling, to interpret coefficients 
# remember pigment_adp and eye_adp are still in original units  
sds 

# summaries
summary(m_tweedie_26)
summary(m_tweedie_45)
summary(m_tweedie_70)

# get response curves for covariates
# these show fitted response + confidence band, on the original (scaled) response scale
cowplot::plot_grid(
  plot( ggeffects::ggpredict( m_tweedie_26, "Body_length_avg")),
  plot( ggeffects::ggpredict( m_tweedie_26, "Femur_elongation")),
  plot( ggeffects::ggpredict( m_tweedie_26, "pigment_adp")),
  plot( ggeffects::ggpredict( m_tweedie_26, "eye_adp"))
)
cowplot::plot_grid(
  plot( ggeffects::ggpredict( m_tweedie_45, "Body_length_avg")),
  plot( ggeffects::ggpredict( m_tweedie_45, "Femur_elongation")),
  plot( ggeffects::ggpredict( m_tweedie_45, "pigment_adp")),
  plot( ggeffects::ggpredict( m_tweedie_45, "eye_adp"))
)
cowplot::plot_grid(
  plot( ggeffects::ggpredict( m_tweedie_70, "Body_length_avg")),
  plot( ggeffects::ggpredict( m_tweedie_70, "Femur_elongation")),
  plot( ggeffects::ggpredict( m_tweedie_70, "pigment_adp")),
  plot( ggeffects::ggpredict( m_tweedie_70, "eye_adp"))
)

# get coefficients with confidence intervals, in response scale (exp())
# for "pigment_adp" and "eye_adp": one ordinal factor increase...
# for "Body_length_avg" and "Femur_elongation": one SD increase...
# ...in x is associated with somewhere between a ... and a ... multiplicative decrease/increase in expected retained
coefs_26 <- as.data.frame( exp( confint(m_tweedie_26, parm = "beta_"))) ; coefs_26
coefs_45 <- as.data.frame( exp( confint(m_tweedie_45, parm = "beta_"))) ; coefs_45
coefs_70 <- as.data.frame( exp( confint(m_tweedie_70, parm = "beta_"))) ; coefs_70

# make rownames a column 
coefs_26 <- rownames_to_column(coefs_26, var = "term")
coefs_45 <- rownames_to_column(coefs_45, var = "term")
coefs_70 <- rownames_to_column(coefs_70, var = "term")

coefs <- cbind( 
  c( 
    rep("SSP1-2.6 2080", nrow(coefs_26)), 
    rep("SSP2-4.5 2080", nrow(coefs_45)), 
    rep("SSP3-7.0 2080", nrow(coefs_70))
    ),
  rbind( coefs_26, coefs_45, coefs_70)
  )

# change column names
names(coefs) <- c("scenario", "term", "lower","upper","estimate")

# remove intercept 
coefs <- subset(coefs, term != "(Intercept)")

# change term names
coefs$term <- rep( c("Body length\n(SD=3.01)", "Leg elongation\n(SD=0.40)", "Depigmentation\n(unscaled)", "Eye regression\n(unscaled)"), 3)

# create directory
dir.create(here::here("4.traits","output"), recursive = TRUE, showWarnings = FALSE)

# export resulting dataframe
save( coefs, file = here::here("4.traits","output","coefs.Rdata"))



# LOSS OF DIVERSITY -----------------------------------------------------------------------------------------------------------------------------------

## CREATING COMMUNITY RASTERS -----------------------------------------------------------------------------------------------------------------------------------

# get suitability projections file names
files <- list.files(here::here("3.range_bagging","output","suitability_rasters"), full.names = TRUE) 

# filter by Boyce, eliminating rows whose stats are based on poor-performing models
files <- files[ keep ]

# quick sanity check
stopifnot(identical(basename(files), paste0("suitability_", stats_sub$species, ".tif")))

# fetch species' suitability prediction rasters
suit_collection <- sprc( files )

# remove "suitability_" and ocasional ".tif" substring from stack names
names(suit_collection) <- sapply( names(suit_collection), function(name) gsub("suitability_", "", name) )
names(suit_collection) <- sapply( names(suit_collection), function(name) gsub("\\.tif", "", name) )

# fetch resolution of raster with finest resolution
# will usually be 1000 by 1000
resolutions <- sapply(suit_collection, function(r) prod(res(r)))
res_target <- res(suit_collection[which.min(resolutions)])

# make empty template grid (biggest extent)
r0 <- rast(ext = ext(suit_collection), 
           resolution = res_target, 
           crs =  crs(suit_collection[1]))

# resample all raster stacks to same template grid
suit_collection_resamp <- lapply( suit_collection, function(r) {
  resample(r, r0, method="bilinear")
})

# get dates/scenarios
dates_scenarios <- names( suit_collection[1] )

# empty list to store output stacks
rasters <- list()

# re-organize rasters to get one raster for each date/scenario, with each species presence/absence prediction as a layer
for (i in dates_scenarios){
  
  # get template grid
  r <- r0
  
  # stack all layers from different species, for same data/scenario
  for (j in 1:length(suit_collection_resamp)){
    r <- c(r, suit_collection_resamp [[j]] [[i]] )
  }
  
  # rename raster layers with species names
  names(r) <- names(suit_collection_resamp)
  
  # output raster to list
  rasters[[length(rasters)+1]] <- r
}

# binarize suitability prediction rasters
rasters_bin <- lapply( rasters, function(r) r > stats_sub$threshold )

# create no-dispersion version of community rasters
# condition: only allow TRUE if TRUE (1) in historical and TRUE (1) in future
rasters_bin_nodisp <- lapply(rasters_bin, function(r){
  r <- rasters_bin[[1]] & r
})

# convert T/F to 1/0
rasters_bin <- lapply( rasters_bin, function(r) r*1 )
rasters_bin_nodisp <- lapply( rasters_bin_nodisp, function(r) r*1 )

# create directory
dir.create(here::here("4.traits","output","community_rasters"), recursive = TRUE, showWarnings = FALSE)

# save community rasters
raster_paths <- here::here("4.traits","output","community_rasters", paste0("communities_", dates_scenarios, ".tif"))
save(raster_paths, file= here::here("4.traits","output","community_rasters","raster_paths.Rdata"))
for( i in seq_along( rasters_bin )){
  writeRaster( rasters_bin[[i]], raster_paths[i], overwrite=TRUE )
}

# save no-dispersion community rasters
raster_paths_nodisp <- here::here( "4.traits","output","community_rasters", paste0("communities_nodisp_", dates_scenarios, ".tif"))
save(raster_paths_nodisp, file=here::here("4.traits","output","community_rasters","raster_paths_nodisp.Rdata"))
for( i in seq_along( rasters_bin_nodisp )){
  writeRaster( rasters_bin_nodisp[[i]], raster_paths_nodisp[i], overwrite=TRUE )
}



## LOADING COMMUNITY RASTERS -----------------------------------------------------------------------------------------------------------------------------------

load(here::here("4.traits","output","community_rasters","raster_paths.Rdata"))
rasters_bin <- lapply( raster_paths, rast )

load(here::here("4.traits","output","community_rasters","raster_paths_nodisp.Rdata"))
rasters_bin_nodisp <- lapply( raster_paths_nodisp, rast )



## TAXONOMIC DIVERSITY -----------------------------------------------------------------------------------------------------------------------------------

# compute taxonomic alpha diversity / species richness in space, for each date/scenario
richness <- lapply( rasters_bin_nodisp, function(r) terra::app(r, sum, na.rm=T) )

# convert list to stack
richness <- do.call(c, richness)
names(richness) <- dates_scenarios

# richness means by data/scenario
lapply( richness, global, na.rm=T)

# plot richness
plot( richness )

# save richness raster
writeRaster( richness, here::here("4.traits","output","richness.tif"), overwrite=TRUE )

# get richness loss by subtracting each future date/scenario to historical 
richness_loss <- richness[[1]] - richness[[-1]]
names(richness_loss) <- dates_scenarios[-1]

# richness loss means by data/scenario
lapply( richness_loss, global, na.rm=T)

# plot richness loss 
plot(richness_loss)
# this map for reference scenario is completely redundant with the "impact" map, although calculating methods differ

# save richness loss raster
writeRaster( richness_loss, here::here("4.traits","output","richness_loss.tif"), overwrite=TRUE )



## FUNCIONAL DIVERSITY -----------------------------------------------------------------------------------------------------------------------------------

### FURTHER PROCESS TRAITS DATA -----------------------------------------------------------------------------------------------------------------------------------

# simplify initial dataframe
traits2 <- traits_all[ c(4, 14:55)]

# subset to species which were in vulnerability analysis and made it through Boyce filter
traits2 <- traits2[ traits2$Genus_species %in% stats_sub$species ,]

# plot continuous traits to check distributions and colinearity
GGally::ggpairs( as.data.frame( traits2[ , sapply( traits2, is.numeric )]))

# checking for outliers with dotplots
par(mfrow=c(5,5)) #TUNE
for (col_name in names( traits2[ , sapply( traits2, is.numeric )] )) {
  dotchart(
    traits2[[col_name]], 
    main = col_name
  )
}

# plotting factor traits and checking for unbalance between classes
facs <- traits2[, sapply(traits2, is.factor)]
par(mfrow = c(5, 5), mar = c(2, 2, 2, 2)) #TUNE
for (i in seq_along(facs)) {
  plot(facs[[i]],  main = names(facs)[i])
}

# select traits to reduce matrix dimensionality and colinearity
traits2 <- traits2[ c("Genus_species", 
                      "Guild", 
                      "Verticality",          
                      "Pigment",
                      "Eyeless",
                      "Eyes_regression",
                      "Body_length_avg", 
                      "Prosoma_length_avg",
                      "Femur_I_length_avg", 
                      "Femur_elongation",
                      "Tibia_I_length_avg",
                      "Capture_web","No_web","Tube_web","Sheet_web","Space_web","Orb_web",               
                      "Active_hunter" )]

# convert species column into row names
rownames(traits2) <- traits2[,1]
traits2 <- traits2[,-1]



### GET FUNCIONAL DISTANCES TREE -----------------------------------------------------------------------------------------------------------------------------------

# create distance matrix with Gower (gower distance is better for mixed data types / data with factors)
dists = BAT::gower( traits2, st = "range")

# build functional tree/dendrogram from distances
func_tree = BAT::tree.build( dists, func = "best" )
# "The best tree is given by bionj with a quality of 0.96"

# plot functional tree/dendrogram
plot(func_tree, cex = 0.5)



### MAP FUNCTIONAL RICHNESS -----------------------------------------------------------------------------------------------------------------------------------

# compute functional diversity in space, for each date/scenario

# standard workflow, but fails for huge rasters
##fd_rasters <- lapply( rasters_bin_nodisp, function(r){
##  
##  # convert raster to cell × species matrix
##  m <- as.matrix(r)
##  
##  # remove cells that are all NA
##  keep2 <- rowSums(!is.na(m)) > 0
##  m <- m[keep2, ]
##  
##  # replace NA with 0 (assuming NA = absence) 
##  m[is.na(m)] <- 0
##  
##  # convert into sparse matrix
##  comm <- Matrix(m, sparse = TRUE)
##  
##  # calculate Faith's PD for all cells
##  fd <- phyloregion::PD(comm, func_tree)
##  
##  # make output raster
##  fd_raster <- r[[1]]
##  values(fd_raster) <- NA
##  values(fd_raster)[keep2] <- fd
##  names(fd_raster) <- "FD"
##  
##  # output
##  fd_raster
##})

# workflow for huge rasters, processes them by blocks
func_diversity_function <- function(r, tree){
  
  # get blocks
  bs <- blocks(r)
  
  # declare empty vector to store results
  fd_vals <- numeric( ncell(r) )
  
  # process by block
  for(i in seq_len( bs$n )) {
    
    # raster to dense matrix
    m <- values(
      r,
      row   = bs$row[i],
      nrows = bs$nrows[i],
      mat   = TRUE
    )
    
    # replace NA with 0 (assuming NA = absence) 
    m[is.na(m)] <- 0
    
    # convert into sparse matrix
    comm <- Matrix::Matrix(m, sparse = TRUE)
    
    # calculate Faith's PD for all cells
    fd <- phyloregion::PD(comm, tree)
    
    # get cell indexes respective to this block
    start <- (bs$row[i] - 1) * ncol(r) + 1
    end   <- start + nrow(m) - 1
    
    # output values in correct vector location
    fd_vals[start:end] <- fd
  }
  
  # make output raster
  fd_rast <- r[[1]]
  values(fd_rast) <- NA
  values(fd_rast) <- fd_vals
  
  # mask output raster
  fd_rast <- crop(fd_rast, richness[[1]], mask=TRUE) 
  
  #output
  fd_rast
}

func_diversity <- lapply( rasters_bin_nodisp, func_diversity_function, tree=func_tree)

# convert list to stack
func_diversity <- do.call(c, func_diversity)
names(func_diversity) <- dates_scenarios
varnames(func_diversity) <- ""

# fd means by data/scenario
lapply( func_diversity, global, na.rm=T)

# plot functional diversity
plot( func_diversity )

# save functional diversity raster
writeRaster( func_diversity, here::here("4.traits","output","func_diversity.tif"), overwrite=TRUE )

# get functional diversity loss by subtracting each future date/scenario with historical 
func_diversity_loss <- func_diversity[[1]] - func_diversity[[-1]]
names(func_diversity_loss) <- dates_scenarios[-1]

# fd loss means by data/scenario
lapply( func_diversity_loss, global, na.rm=T)

# plot functional diversity loss 
plot( func_diversity_loss )

# save functional diversity loss raster
writeRaster( func_diversity_loss, here::here("4.traits","output","func_diversity_loss.tif"), overwrite=TRUE )

# get percentage functional diversity loss
per_func_diversity_loss <- func_diversity_loss / (func_diversity[[1]]+0.00000001) * 100

# save percentage functional diversity loss raster
writeRaster( per_func_diversity_loss, here::here("4.traits","output","per_func_diversity_loss.tif"), overwrite=TRUE )

# portion of "europe" with more than 25% and 50% functional diversity loss
over_25 <- lapply( per_func_diversity_loss, function(r) global(r > 25, "mean", na.rm = TRUE)) ; over_25 #RESULT
over_50 <- lapply( per_func_diversity_loss, function(r) global(r > 50, "mean", na.rm = TRUE)) ; over_50 #RESULT



### DISCERN FUNCTIONAL DIVERSITY LOSS FROM RICHNESS LOSS -----------------------------------------------------------------------------------------------------------------------------------
# create a null distribution of functional diversity loss values
# which sites have higher or lower functional diversity loss than expected by chance, based on richness loss?

# creating template rasters
random_FDlosses_26_2080 = func_diversity_loss[[2]]
random_FDlosses_45_2080 = func_diversity_loss[[4]]
random_FDlosses_70_2080 = func_diversity_loss[[6]]

# iterating 
for (i in 1:100){
  cat(paste("Run", i, "\n"))
  
  # getting random functional distance tree by shuffling species names
  random_tree <- func_tree
  random_tree$tip.label = sample(func_tree$tip.label)
  
  # computing funcional diversity
  # rasters_bin_nodisp[c(1,3,5,7)] means subsetting only historical and 2080 scenarios
  random_func_diversity = lapply( rasters_bin_nodisp[c(1,3,5,7)], func_diversity_function, tree=random_tree )
  
  # convert list to stack
  random_func_diversity <- do.call(c, random_func_diversity)

  # get functional diversity loss by subtracting each future date/scenario with historical 
  random_func_diversity_loss <-  random_func_diversity[[1]] - random_func_diversity[[-1]]

  random_FDlosses_26_2080 <- c(random_FDlosses_26_2080, random_func_diversity_loss[[1]])
  random_FDlosses_45_2080 <- c(random_FDlosses_45_2080, random_func_diversity_loss[[2]])
  random_FDlosses_70_2080 <- c(random_FDlosses_70_2080, random_func_diversity_loss[[3]])
}

# removing template from each grid
random_FDlosses_26_2080 <- random_FDlosses_26_2080[[-1]]
random_FDlosses_45_2080 <- random_FDlosses_45_2080[[-1]]
random_FDlosses_70_2080 <- random_FDlosses_70_2080[[-1]]

# save 
writeRaster(random_FDlosses_26_2080, here::here("4.traits","output","random_FDlosses_26_2080.tif"), overwrite=TRUE )
writeRaster(random_FDlosses_45_2080, here::here("4.traits","output","random_FDlosses_45_2080.tif"), overwrite=TRUE )
writeRaster(random_FDlosses_70_2080, here::here("4.traits","output","random_FDlosses_70_2080.tif"), overwrite=TRUE )

# calculate and map SES 26_2080
ses_26_2080 = (func_diversity_loss[[2]] - mean( random_FDlosses_26_2080 )) / app(random_FDlosses_26_2080, fun = sd)
ses_26_2080 ; global(ses_26_2080, na.rm=TRUE) ; plot(ses_26_2080)

# calculate and map SES 45_2080
ses_45_2080 = (func_diversity_loss[[4]] - mean( random_FDlosses_45_2080 )) / app(random_FDlosses_45_2080, fun = sd)
ses_45_2080 ; global(ses_45_2080, na.rm=TRUE) ; plot(ses_45_2080)

# calculate and map SES 70_2080
ses_70_2080 = (func_diversity_loss[[6]] - mean( random_FDlosses_70_2080 )) / app(random_FDlosses_70_2080, fun = sd)
ses_70_2080 ; global(ses_70_2080, na.rm=TRUE) ; plot(ses_70_2080)

# high values denote higher functional diversity loss than expected based on richness loss
# low values denote lower functional diversity loss than expected based on richness loss
# places with sd = 0 come out NA, make them gray

# save
writeRaster(ses_26_2080, here::here("4.traits","output","ses_26_2080.tif"), overwrite=TRUE )
writeRaster(ses_45_2080, here::here("4.traits","output","ses_45_2080.tif"), overwrite=TRUE )
writeRaster(ses_70_2080, here::here("4.traits","output","ses_70_2080.tif"), overwrite=TRUE )

