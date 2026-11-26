
# These are the scripts to see the range bagging results and create all figures.

library( here )
library( readxl )
library( ggplot2 )
library(cowplot)
library( dplyr )
library( tidyr )
library( forcats )
library( ggspatial ) # plotting maps
library( terra )
library( raster )
library(ggimage)
library(investr)
library(ggpubr)
library(ggbeeswarm)
library(gghalves)
library(plyr)

set.seed(46758)

# create directory
dir.create(here::here("5.results_and_figures","figures"), recursive = TRUE, showWarnings = FALSE)



# IMPORTING DATA -----------------------------------------------------------------------------------------------------------------------------------

# load output object
load(here::here("3.range_bagging","output","stats.RData"))

# convert all columns except species to numeric
stats[,-1] <- lapply(stats[,-1], as.numeric)

# filter by Boyce, eliminating rows whose stats are based on poor-performing models
nrow(stats) ; stats <- subset( stats, boyce > 0.5) ; nrow(stats) #TUNE boyce threshold 

# import functional traits database
traits <- read_excel( here::here("4.traits","Database_Mammola_et_al_2022_Figshare_V3.xlsx"))

# uniformizing species name structure
traits$Genus_species <- gsub( " ", "_", traits$Genus_species)
traits$Genus_species <- gsub( "\\.", "", traits$Genus_species)
traits[ traits$Genus_species == "Gen_sp_1", "Genus_species" ] <- "Dysderidae_Gen."
traits[ traits$Genus_species == "Gen_12_sp_1", "Genus_species" ] <- "Linyphiidae_Gen."



# MAIN RESULTS TABLE -----------------------------------------------------------------------------------------------------------------------------------

# get stat means and sd in readable format
data.frame( 
  Mean = format( colMeans(stats[,-1], na.rm = TRUE), scientific = FALSE ),
  SD = format( sapply( stats[,-1], sd, na.rm = TRUE), scientific = FALSE )
) #RESULT

# amount of species with total site loss
count_total_site_loss <- apply( 
  stats[, grepl( "p_site_loss", colnames(stats))] ,
  2, function(x) table(x)[names(table(x))==1]
  )
count_total_site_loss #RESULT

# percentage of species of total site loss
count_total_site_loss / nrow(stats) #RESULT



# RAREFACTION CURVES FIGURE  -----------------------------------------------------------------------------------------------------------------------------------

## ALL-SPECIES RAREFACTION CURVES PLOT -----------------------------------------------------------------------------------------------------------------------------------

# load output object
load(here::here("3.range_bagging","output","rarefaction_plot_data_d.RData"))

# combine all iterations
all_preds <- bind_rows(rarefaction_plot_data_d)

# create empty column to store class
all_preds$class <- character(length= nrow(all_preds))

# add ecological classification of each species to dataframe
for (i in 1:nrow(all_preds)) {
  shortname <- all_preds[i,3]
  all_preds[i,4] <- traits[ traits$Genus_species == shortname, "Ecological_classification"]
}

# plot
rarefaction_plot <- ggplot( 
  all_preds,
  aes(
    x = n_locations,
    y = mean_vertices_dist,
    color = class,
    group = species 
  )
) + 
  geom_line( linewidth = 0.3) + 
  xlab("Number of sampled locations") +
  ylab("Captured niche") + 
  scale_x_continuous(trans='log10') +
  theme_classic()+
  
  scale_color_manual(
    name = "Ecological\nclassification", # legend title
    values = c("azure4", "burlywood3")
  )+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.2) ,
    legend.text = element_text(size = 11)
  )


## EXAMPLE OF SINGLE RAREFACTION CURVE  -----------------------------------------------------------------------------------------------------------------------------------

# import example rarefaction data and fit from "Amaurobius ferox" - troglophile
load( here::here("3.range_bagging","output","example_rarefaction.Rdata") )
load( here::here("3.range_bagging","output","example_weibull_fit.Rdata") )

# get values
real_niche_width_d <- summary(exp_fit)$parameters[1,1]
sampled_niche_width_d <- max(rarified_d$mean_vertices_dist)
inflate_d <-  real_niche_width_d / sampled_niche_width_d

# x values for plot
xs <- seq(3.5, 140, length.out = 100)

# get predicted values from fit
preds <- investr::predFit(
  exp_fit,
  newdata = data.frame(n_locations = xs),
  interval = "confidence"
)

# restructure
pred_df <- data.frame(
  n_locations = xs,
  fit = preds[, "fit"],
  lwr = preds[, "lwr"],
  upr = preds[, "upr"]
)

ex1 <- ggplot() +
  # confidence band
  geom_ribbon(
    data = pred_df,
    aes(x = n_locations, ymin = lwr, ymax = upr),
    fill = "burlywood2",
    alpha = 0.3
  ) +
  # fitted line
  geom_line(
    data = pred_df,
    aes(x = n_locations, y = fit),
    color = "burlywood3",
    linewidth = 1
  ) +
  # observed data
  geom_point(
    data = rarified_d,
    aes(x = n_locations, y = mean_vertices_dist),
    size = 0.5
  ) +
  theme_classic() +
  labs(
    x = "Number of sampled locations",
    y = NULL
  ) +
  annotate(
    "text",
    x = 65, y = 1.3,
    label = paste(
      paste0("Sampled niche: ", round( sampled_niche_width_d, 2 )),
      paste0("Real niche: ", round( real_niche_width_d, 2 ), "±", round( summary(exp_fit)$parameters[1,2], 2)),
      paste0("Expansion factor: ", round( inflate_d, 2 )), 
      sep = "\n"),
    hjust = 0,
    vjust = 1,
    size = 4
  ) 


# another from "Folkia inermis" - troglobiont
load( here::here("3.range_bagging","output","example_rarefaction2.Rdata") )
load( here::here("3.range_bagging","output","example_weibull_fit2.Rdata") )

# get values
real_niche_width_d2 <- summary(exp_fit)$parameters[1,1]
sampled_niche_width_d2 <- max(rarified_d$mean_vertices_dist)
inflate_d2 <-  real_niche_width_d2 / sampled_niche_width_d2

# x values for plot
xs2 <- seq(3.5, 14, length.out = 100)

# get predicted values from fit
preds2 <- investr::predFit(
  exp_fit,
  newdata = data.frame(n_locations = xs2),
  interval = "confidence"
)

# restructure
pred_df2 <- data.frame(
  n_locations = xs2,
  fit = preds2[, "fit"],
  lwr = preds2[, "lwr"],
  upr = preds2[, "upr"]
)

ex2 <- ggplot() +
  # confidence band
  geom_ribbon(
    data = pred_df2,
    aes(x = n_locations, ymin = lwr, ymax = upr),
    fill = "azure3",
    alpha = 0.3
  ) +
  # fitted line
  geom_line(
    data = pred_df2,
    aes(x = n_locations, y = fit),
    color = "azure4",
    linewidth = 1
  ) +
  # observed data
  geom_point(
    data = rarified_d,
    aes(x = n_locations, y = mean_vertices_dist),
    size = 0.8
  ) +
  theme_classic() +
  labs(
    x = "Number of sampled locations",
    y = NULL
  ) +
  annotate(
    "text",
    x = 9, y = 0.228,
    label = paste(
      paste0("Sampled niche: ", round( sampled_niche_width_d2, 2 )),
      paste0("Real niche: ", round( real_niche_width_d2, 2 ), "±", round( summary(exp_fit)$parameters[1,2], 2)),
      paste0("Expansion factor: ", round( inflate_d2, 2 )), 
      sep = "\n"),
    hjust = 0,
    vjust = 1,
    size = 4
    )



## JOIN ALL AND SINGLE RAREFACTION PLOTS TO MAKE FIGURE -----------------------------------------------------------------------------------------------------------------------------------

# define plot grid
fig <- plot_grid( 
  NULL, rarefaction_plot, NULL, ex1, NULL, ex2, 
  nrow = 1, 
  rel_widths = c(0.15, 1, 0.15, 1, 0.15, 1))

# add remaining elements
fig <- ggdraw(fig) +
  
  draw_plot_label(
    c("a)", "b)", "c)"),
    x = c(0.04, 0.36, 0.703),
    y = c(0.90, 0.90, 0.90),
    hjust = 0, vjust = 0)

# export
ggsave( plot = fig, 
        filename = here::here("5.results_and_figures","figures","rarefaction.png"), path = NULL,
        width = 16, height = 4,
        bg = "white",
        dpi=800)



# HABITAT SUITABILITY, SITE LOSS AND RANGE CHANGE FIGURE -----------------------------------------------------------------------------------------------------------------------------------

# refine dataset
data <- stats %>%
  dplyr::select(
    species,
    mean_suit_historical,
    mean_suit_future_26_2080,
    mean_suit_future_45_2080,
    mean_suit_future_70_2080
  ) %>%
  tidyr::pivot_longer(
    cols = -species,
    names_to = "text",
    values_to = "value"
  ) %>%
  mutate(text = fct_reorder(text, value))

# create empty column to store class
data$class <- character(length= nrow(data))

# add ecological classification of each species to dataframe
for (i in 1:nrow(data)) {
  shortname <- data[i,1]
  data[i,4] <- traits[ traits$Genus_species == as.character(shortname), "Ecological_classification"]
}

# summary statistics to plot mean and confidence interval
result_tot <- data %>%
  dplyr::group_by(text, class) %>%
  dplyr::summarise(
    ES = mean(value, na.rm = TRUE), # effect size / mean
    SD = sd(value, na.rm = TRUE), # standard deviation
    n = dplyr::n(), 
    SE = SD / sqrt(n), # standard error
    L = ES - 1.96 * SE, # lower limit of 95% confidence interval
    U = ES + 1.96 * SE, # upper limit of 95% confidence interval
    .groups = "drop"
  )

# get t-Test results
# test differences in distributions/mean of suitability with paired Student's t-Test
# will suitability decrease in the future (historical vs futures)?
future_suits <- list( stats$mean_suit_future_26_2080,
                      stats$mean_suit_future_45_2080,
                      stats$mean_suit_future_70_2080 )
t_results <- lapply( future_suits, t.test, stats$mean_suit_historical, paired = TRUE)

# organize results in dataframe
t_df <- data.frame(
  text = c( "mean_suit_future_26_2080", "mean_suit_future_45_2080", "mean_suit_future_70_2080"),
  t_statistic = sapply(t_results, function(x) unname(x$statistic)),
  p_value = sapply(t_results, function(x) x$p.value)
)

# function to get unicode exponents
to_superscript <- function(x) {
  supers <- c(
    "0"="⁰", "1"="¹", "2"="²", "3"="³", "4"="⁴",
    "5"="⁵", "6"="⁶", "7"="⁷", "8"="⁸", "9"="⁹",
    "-"="⁻"
  )
  
  s <- format(x, scientific = TRUE)
  parts <- strsplit(s, "e")[[1]]
  mantissa <- sprintf("%.2f", as.numeric(parts[1]))
  exponent <- sub("^\\+", "", parts[2])
  exponent <- sub("^(-?)0+", "\\1", exponent)
  
  exp_sup <- paste0(supers[strsplit(exponent, "")[[1]]], collapse = "")
  paste0(mantissa, "×10", exp_sup)
}

# construct labels that concatenate t_statistic + p_value
t_df$label <- sprintf(
  "t = %.2f\np = %s",
  t_df$t_statistic,
  sapply( t_df$p_value, to_superscript )
)

# declare plot
suit_plot <- ggplot(data, aes(
  x = forcats::fct_rev(text),
  y = value,
  fill = class
)) +
  geom_half_violin(
    data = subset(data, class == "Troglobiont"),
    side = "l",
    color = NA,
    width = 0.4,
    trim = FALSE
  ) +
  geom_half_violin( 
    data = subset(data, class == "Troglophile"),
    side = "r",
    color = NA,
    width = 0.4,
    trim = FALSE
  ) +

  scale_fill_manual(
    values = c(
      Troglobiont = "azure3",
      Troglophile = "burlywood3"
    )
  ) +
  
  # data points
  geom_jitter(
    data = subset(data, class == "Troglobiont"),
    aes(x = as.numeric(forcats::fct_rev(text)) - 0.3, y = value),
    color = "azure3",
    shape = 16,
    size = 1,
    width = 0.05,
    alpha = 1
  ) +
  geom_jitter(
    data = subset(data, class == "Troglophile"),
    aes(x = as.numeric(forcats::fct_rev(text)) + 0.3, y = value),
    color = "burlywood3",
    shape = 16,
    size = 1,
    width = 0.05,
    alpha = 1
  ) +
  
  # mean and confidence interval
  geom_pointrange(
    data = subset(result_tot, class == "Troglobiont"),
    aes(x = forcats::fct_rev(text), y = ES, ymin = L, ymax = U),
    size = 0.7,
    color = "black",
    position = position_nudge(x = -0.3),
    show.legend = FALSE
  ) +
  geom_pointrange(
    data = subset(result_tot, class == "Troglophile"),
    aes(x = forcats::fct_rev(text), y = ES, ymin = L, ymax = U),
    size = 0.7,
    color = "black",
    position = position_nudge(x = 0.3),
    show.legend = FALSE
  ) +

  theme_classic() +
  
  theme(legend.position = "right", text = element_text(size = 13)) +
  
  labs( x = "Scenario and period", y = "Mean habitat suitability inside range" ) +
  
  scale_fill_manual(
    labels = c(
      Troglobiont = "Troglobionts",
      Troglophile = "Troglophiles"
    ),
    values = c(
      Troglobiont = "azure3",
      Troglophile = "burlywood3"
    ),
    name = "Ecological\nclassification"
  ) +
  
  scale_x_discrete(
    labels = c(
    "mean_suit_historical" = "Historical",
    "mean_suit_future_26_2080" = "SSP1-2.6 2080",
    "mean_suit_future_45_2080" = "SSP2-4.5 2080",
    "mean_suit_future_70_2080" = "SSP3-7.0 2080"
    )
  ) +
  
  # t_test results
  geom_text(
    data = t_df,
    aes(
      x = forcats::fct_rev(text),
      y = 1.2,
      label = label
    ),
    inherit.aes = FALSE
  )



## SITE LOSS PLOT -----------------------------------------------------------------------------------------------------------------------------------

# refine dataset
data <- stats %>%
  dplyr::select(
    species,
    p_site_loss_future_26_2080,
    p_site_loss_future_45_2080,
    p_site_loss_future_70_2080
  ) %>%
  tidyr::pivot_longer(
    cols = -species,
    names_to = "text",
    values_to = "value"
  ) 

# create empty column to store class
data$class <- character(length= nrow(data))

# add ecological classification of each species to dataframe
for (i in 1:nrow(data)) {
  shortname <- data[i,1]
  data[i,4] <- traits[ traits$Genus_species == as.character(shortname), "Ecological_classification"]
}

# plot
p_site <- ggplot() +
  
  geom_quasirandom( 
    data = subset(data, class == "Troglophile"),
    aes(x = as.numeric(as.factor(text)) - 0.15, y = value),
    alpha = 1, 
    width = 0.04, 
    colour = "burlywood3",
    cex = 1, 
    show.legend = FALSE, 
    bandwidth = 2) +
  geom_quasirandom( 
    data = subset(data, class == "Troglobiont"),
    aes(x = as.numeric(as.factor(text)) + 0.15, y = value),
    alpha = 1, 
    width = 0.04, 
    colour = "azure3",
    cex = 1, 
    show.legend = FALSE, 
    bandwidth = 2) +
  
  geom_bar(
    data = subset(data, class == "Troglophile"),
    aes(x = text , y = value),
    stat="mean", 
    width=0.2, 
    fill = "burlywood3", 
    alpha = 0.6, 
    position = position_nudge(x = -0.15)) +
  geom_bar(
    data = subset(data, class == "Troglobiont"),
    aes(x = text , y = value),
    stat="mean", 
    width=0.2, 
    fill = "azure3", 
    alpha = 0.6, 
    position = position_nudge(x = +0.15)) +
  
  labs( x = "Scenario and period", y = "Fraction of current sites predicted to be lost") +
  
  theme_classic()+
  
  scale_x_discrete(
    labels = c(
      "p_site_loss_future_26_2080" = "SSP1-2.6 2080",
      "p_site_loss_future_45_2080" = "SSP2-4.5 2080",
      "p_site_loss_future_70_2080" = "SSP3-7.0 2080"
    )
  )



## RANGE CHANGE STATISTICS TERNARY PLOT --------------------------------------------------------------------------------------------------------

scen <- c("26" = "SSP1-2.6 2080",
          "45" = "SSP2-4.5 2080",
          "70" = "SSP3-7.0 2080")

tern_data <- do.call(rbind, lapply(names(scen), function(code) {
  data.frame(
    species  = stats$species,
    scenario = scen[[code]],
    gain     = stats[[paste0("p_gained_future_", code, "_2080")]],
    loss     = -stats[[paste0("p_lost_future_",  code, "_2080")]]
  )
}))

tern_data <- tern_data %>%
  tidyr::drop_na(gain, loss) %>%
  dplyr::mutate(
    stable = 1 - loss,
    total  = gain + loss + stable,
    Gain   = gain   / total,
    Loss   = loss   / total,
    Stable = stable / total,
    net    = gain - loss
  )

# add ecological classification 
tern_data$class <- traits$Ecological_classification[match(tern_data$species, traits$Genus_species)]
tern_data <- dplyr::filter(tern_data, !is.na(class))

tern_data$scenario <- factor(tern_data$scenario, levels = scen)
tern_data$class    <- factor(tern_data$class, levels = c("Troglophile", "Troglobiont"))

# mean position per SSP and ecological class (centroids)
tern_centroids <- tern_data %>%
  dplyr::group_by(scenario, class) %>%
  dplyr::summarise(
    Gain   = mean(Gain),
    Loss   = mean(Loss),
    Stable = mean(Stable),
    .groups = "drop"
  )


library(ggtern)

p_tern <- ggtern(tern_data, aes(x = Loss, y = Gain, z = Stable)) +
  
  geom_point(aes(colour = net, shape = class), size = 1.8, alpha = 0.85) +
  
  scale_colour_gradient2(
    name = "Net range\nchange",
    low = "firebrick", mid = "grey85", high = "steelblue",
    midpoint = 0, limits = c(-1, 1)
  ) +
  scale_shape_manual(
    name = "Ecological\nclassification",
    values = c(Troglophile = 16, Troglobiont = 17),
    labels = c(Troglophile = "Troglophiles", Troglobiont = "Troglobionts")
  ) +
  guides(shape = guide_legend(override.aes = list(colour = "grey30", size = 3))) +
  
  # centroids: same symbol as the group, filled black
  geom_point(data = subset(tern_centroids, class == "Troglophile"),
             shape = 21, fill = "black", colour = "white", size = 4.5, stroke = 0.4) +
  geom_point(data = subset(tern_centroids, class == "Troglobiont"),
             shape = 24, fill = "black", colour = "white", size = 5, stroke = 0.4) +
  
  facet_wrap(~ scenario, nrow = 1, strip.position = "bottom") +
  labs(x = "Loss", y = "Gain", z = "Stable") +
  
  theme_bw() +
  theme_showarrows() +
  theme_nomask() +
  theme(
    strip.background = element_blank(),     
    strip.text       = element_text(size = 12, face = "plain"),
    strip.placement = "outside",
    panel.spacing.x  = unit(1.5, "cm"), 
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 7)
  )


## MAKE FIGURE -----------------------------------------------------------------------------------------------------------------------------------

# extract legend from suitability plot
legend <- get_legend( suit_plot )

# ggtern plots must be converted to a grob before cowplot can place them
p_tern_grob <- ggtern::ggplotGrob( p_tern )

# top row: suitability + site loss
top_row <- plot_grid(
  NULL, suit_plot + theme(legend.position = "none"), NULL, p_site, NULL,
  nrow = 1, rel_widths = c(0.1, 1.0, 0.2, 0.6, 0.3)
)

# bottom row: ternary plot
bottom_row <- plot_grid(
  NULL, p_tern_grob, NULL,
  nrow = 1, rel_widths = c(0.03, 1, 0.03)
)

# stack rows
fig <- plot_grid(
  top_row, bottom_row,
  nrow = 2, rel_heights = c(1, 1.1)
) +
  
  draw_plot_label(
    c("a)", "b)", "c)"),
    x = c(0.02, 0.55, 0.02),
    y = c(0.93, 0.93, 0.50),
    hjust = 0, vjust = 0
  ) +
  
  # add legend
  draw_plot(legend, x = 0.895, y = 0.6, width = 0.08, height = 0.3)

# export
ggsave( plot = fig,
        filename = here::here("5.results_and_figures","figures","suitability+site_loss+ternary.png"),
        width = 16, height = 10,
        bg = "white",
        dpi = 600)


# IMPACT AND MESS MAPS FIGURE -----------------------------------------------------------------------------------------------------------------------------------

## EUROPE IMPACT MAP -----------------------------------------------------------------------------------------------------------------------------------

# load impact raster
impact_sum <- rast(here::here("3.range_bagging","output","impact_sum.tif"))

# convert spatraster do df, for ggplot
impact_df <- as.data.frame(impact_sum, xy = TRUE, na.rm = TRUE)
colnames(impact_df)[3] <- "impact"

# declare europe plot
i_plot <- ggplot() +
  geom_raster(data = impact_df,
              aes(x = impact_df$x, y = impact_df$y, fill = impact_df$impact)) +
  scale_fill_gradient(
    low = "#e8e8e8",
    high = "red",
    name = "Richness loss"
  ) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none") +
  
  # north
  annotation_north_arrow(
    location = "tr",
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm"),
    pad_x = unit(2.5, "cm"),
    pad_y = unit(1, "cm"),
    style = north_arrow_minimal
  ) +
  
  # map scale 
  annotation_scale( 
    style = "ticks", 
    location = "tr", 
    pad_x = unit(4.5, "cm"), 
    pad_y = unit(1, "cm")
  ) +
  
  # inset rectangle
  annotate(
    "rect",
    xmin =  -700000, xmax = 500000,
    ymin = 1200000, ymax = 1980000,
    colour = "black",
    fill = NA,
    linewidth = 0.4
  )

# extract legend
legend_i <- get_legend( i_plot + theme(legend.position = "right"))



## DINARICS IMPACT INSET MAP -----------------------------------------------------------------------------------------------------------------------------------

# crop raster to dinarics
impact_sum_din <- crop(impact_sum, ext( -700000, 500000, 1200000, 1910000))

# convert spatraster do df, for ggplot
impact_df_din <- as.data.frame(impact_sum_din, xy = TRUE, na.rm = TRUE)
colnames(impact_df_din)[3] <- "impact"

# declare dinarics plot
i_plot_din <- ggplot() +
  geom_raster(data = impact_df_din,
              aes(x = impact_df_din$x, y = impact_df_din$y, fill = impact_df_din$impact)) +
  scale_fill_gradient(
    low = "#e8e8e8",
    high = "red",
    name = "Impact"
  ) +
  coord_equal() +
  theme_void()+
  theme(legend.position = "none") +
  
  # map scale 
  annotation_scale( 
    style = "ticks", 
    location = "tr", 
    pad_x = unit(1, "cm"), 
    pad_y = unit(1, "cm")
  ) +
  
  # inset rectangle
  annotate(
    "rect",
    xmin =  -700000, xmax = 500000,
    ymin = 1200000, ymax = 1910000,
    colour = "black",
    fill = NA,
    linewidth = 0.4
  )



## MESS MAP -----------------------------------------------------------------------------------------------------------------------------------

# load MESS raster
MESS_mean <- rast(here::here("3.range_bagging","output","MESS_mean.tif"))

# get extrapolation by substituting all positive values (similarity) by zero, since they mean no extrapolation
extrapolation <- -ifel(MESS_mean > 0, 0, MESS_mean)

# convert spatraster do df, for ggplot
extrapolation_df <- as.data.frame(extrapolation, xy = TRUE, na.rm = TRUE)
colnames(extrapolation_df)[3] <- "extrapolation"

# declare plot
e_plot <- ggplot() +
  geom_raster(data = extrapolation_df,
              aes(x = extrapolation_df$x, y = extrapolation_df$y, fill = extrapolation_df$extrapolation)) +
  scale_fill_gradient(
    low = "#e8e8e8",
    high = "purple4",
    trans = scales::pseudo_log_trans(sigma = 10),
    name = "Extrapolation\nLog-scale"
  ) +
  coord_equal() +
  theme_void()+
  theme(legend.position = "none") +
  guides(
    fill = guide_colorbar(
      barheight = unit(8, "cm"),
      barwidth = unit(0.7, "cm")
    ))

# extract legend
legend_e <- get_legend( e_plot + theme(legend.position = "right"))



## IMPACT/MESS BIPLOT -----------------------------------------------------------------------------------------------------------------------------------

raster1 <- raster::raster(here::here("3.range_bagging","output","impact_sum.tif"))
raster2 <- -raster::raster(here::here("3.range_bagging","output","MESS_mean.tif")) # note minus sign: gets inverse of MESS, where higher values are higher extrapolation

# classify each raster into 2 or 3 categories
breaks1 <- c(-0.001, 5, 10, 25)
breaks2 <- c(-100, 0, 600)

raster1_class <- cut(raster1[], breaks = breaks1, labels = c(1, 2, 3))  
raster2_class <- cut(raster2[], breaks = breaks2, labels = c(1, 2))  

# combine rasters
combined_vals <- as.numeric(raster1_class) * 10 + as.numeric(raster2_class)
combined_raster <- terra::rast(here::here("3.range_bagging","output","impact_sum.tif"))
combined_raster[] <- combined_vals

# define palette
bivariate_colors <- c(
  "Low_Low" = "#e8e8e8", "Low_Medium" = "#ace4e4",  # Low values of Raster1
  "Medium_Low" = "#dfb0d6", "Medium_Medium" = "#a5add3",  # Medium values of Raster1
  "High_Low" = "#be64ac", "High_Medium" = "#8c62aa"   # High values of Raster1
)

# convert raster to a data frame for plotting 
raster_df <- as.data.frame(combined_raster, xy = TRUE)
raster_df$layer <- as.character(raster_df$sum)

# assign labels to raster values
raster_df$layer <- plyr::revalue(raster_df$layer, c("11" = "Low_Low", "12" = "Low_Medium",
                                              "21" = "Medium_Low", "22" = "Medium_Medium",
                                              "31" = "High_Low", "32" = "High_Medium"))

# assign colors to labels using the palette
raster_df$color <- bivariate_colors[raster_df$layer] 

# plot map  
biplot <- ggplot() +
  geom_raster(
    data = raster_df, 
    aes(x = x, y = y, fill = color)) +
  scale_fill_identity() +
  coord_equal() +
  theme_void()

# create legend dataframe
legend_df <- expand.grid(
  Extrapolation = c("No", "Yes"),  # Categories for Raster 2
  Impact = c("Low", "Medium", "High")  # Categories for Raster 1
)
legend_df$color <- as.vector(bivariate_colors)  # Add the colors to the data frame

# plot legend
legend_bi <- ggplot() +
  geom_tile(data = legend_df, aes(x = Extrapolation, y = Impact, fill = color)) + # raster 2, raster 1
  scale_fill_identity() +
  labs(x = "Extrapolation →", y = "Richness loss →") + # raster 2, raster 1 LABELS
  theme_minimal() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))



## JOIN ALL FOUR MAPS TO MAKE FIGURE -----------------------------------------------------------------------------------------------------------------------------------

# put dinaric inset inside a grid with a scale, to reduce its size in final figure
i_plot_din_reduced <- plot_grid(i_plot_din, scale = 0.7)

# define plot grid
fig <- plot_grid( 
  i_plot, i_plot_din_reduced, 
  e_plot, biplot, 
  nrow = 2 )

# add remaining elements
fig <- ggdraw(fig) +
  
  # legends
  draw_plot(legend_i, x = 0.025, y = 0.55, width = 0.08, height = 0.3) +
  draw_plot(legend_e, x = 0.025, y = 0.15, width = 0.08, height = 0.3) +
  draw_plot(legend_bi, x = 0.54, y = 0.29, width = 0.1, height = 0.14) +
  
  # inset line
  draw_line(
    x = c(0.333, 0.591),
    y = c(0.700, 0.875),
    colour = "black",
    linewidth = 0.4
  ) +
  
  draw_plot_label(
    c("a)", "b)", "c)"),
    x = c(0.04, 0.04, 0.54),
    y = c(0.92, 0.50, 0.50),
    hjust = 0, vjust = 0)

# export
ggsave( plot = fig, 
        filename = here::here("5.results_and_figures","figures","impact+MESS.png"),
        width = 16, height = 12,
        bg = "white",
        dpi=600)



# TRAITS FIGURE

## TRAIT EFFECT SIZE FOREST PLOT DECOMPOSED BY SCENARIO -----------------------------------------------------------------------------------------------------------------------------------

# import coefs data
load(here::here("4.traits","output","coefs.Rdata"))

# changing factor order for plot
coefs$scenario <- factor(
  coefs$scenario,
  levels = c(
    "SSP3-7.0 2080",
    "SSP2-4.5 2080",
    "SSP1-2.6 2080"
  )
)

# plot
forest <- ggplot(coefs,
                 aes(y = term,
                     x = estimate,
                     colour = scenario)) +
  
  geom_vline(xintercept = 1,
             linetype = "dashed",
             colour = "grey50") +
  
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    position = position_dodge(width = 0.5),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  
  scale_colour_manual(
    values = c(
      "SSP1-2.6 2080" = "aquamarine3",
      "SSP2-4.5 2080" = "coral",
      "SSP3-7.0 2080" = "coral3"
    ),
    breaks = c( # legend order
      "SSP1-2.6 2080",
      "SSP2-4.5 2080",
      "SSP3-7.0 2080"
    ),
    name = "Scenario"
  ) +
  
  scale_x_log10() +
  
  labs(
    x = "Effect size: exp(β) ± 95% CI",
    y = NULL,
    colour = ""
  ) +
  theme_classic()

# x=1 means no effect (multiplies prediction by 1)
# negative effects and negative response values mean negative range change (loss)
# positive effects and positive response values mean less-negative range change (less loss, or even gain in few cases)

# fetching silhouettes
sil_df <- data.frame(
  term = c(
    "Leg elongation\n(SD=0.40)",
    "Eye regression\n(unscaled)",
    "Depigmentation\n(unscaled)",
    "Body length\n(SD=3.01)"
  ),
  image = here::here(
    "4.traits",
    "spider_silhouettes",
    c("g56.png", "g48.png", "g36.png", "g54.png")
  ),
  size = c(0.11, 0.10, 0.15, 0.15)
)
icons <- ggplot(sil_df, aes(x = 1, y = term)) +
  geom_image(aes(image = image, size = size)) +
  scale_size_identity() +
  scale_x_continuous(limits = c(0.8, 1.2)) +
  theme_void()



## SES PLOT -----------------------------------------------------------------------------------------------------------------------------------

# load rasters
richness <- rast(here::here("4.traits","output","richness.tif")) # for mask
ses_45_2080 <- rast(here::here("4.traits","output","ses_45_2080.tif"))

# convert spatraster do df, for ggplot
grey_mask <- as.data.frame(richness[[1]], xy = TRUE, na.rm = TRUE) # for mask
ses_df <- as.data.frame(ses_45_2080, xy = TRUE, na.rm = TRUE)
colnames(ses_df)[3] <- "SES"

# declare plot
ses_plot <- ggplot() +
  
  # grey mask
  geom_raster(
    data = grey_mask,
    aes(x, y),
    fill = "#e8e8e8" 
  ) +
  
  # ses
  geom_raster(data = ses_df,
              aes(x = ses_df$x, y = ses_df$y, fill = ses_df$SES)) +
  scale_fill_gradient2(
    low = "#2166AC",     
    mid = "white",     
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    name = "SES\nFunctional richness loss"
  ) +
  coord_equal() +
  theme_void() +
  theme(legend.position = c(0.96,0.65)) +
  
  # north
  annotation_north_arrow(
    location = "tr",
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm"),
    pad_x = unit(2, "cm"),
    pad_y = unit(1, "cm"),
    style = north_arrow_minimal
  ) +
  
  # map scale 
  annotation_scale( 
    style = "ticks", 
    location = "tl", 
    pad_x = unit(6, "cm"), 
    pad_y = unit(2, "cm")
  ) 



## FINAL FIGURE  -----------------------------------------------------------------------------------------------------------------------------------

# combine icons and forest
forest_icons <- plot_grid(
  icons,
  forest,
  nrow = 1,
  rel_widths = c(0.12, 1),
  align = "h",
  axis = "tb"
)

# shorten forest plot
forest_short <- ggdraw() +
  draw_plot(forest_icons, x = 0, y = 0.15, width = 1, height = 0.7)

# define plot grid
fd_fig <- plot_grid(
  
  NULL, forest_short, ses_plot,
  nrow = 1,
  rel_widths = c(0.2, 1, 2)
)+
  
  draw_plot_label(
    c("a)", "b)"),
    x = c(0.02, 0.45),
    y = c(0.93, 0.92),
    hjust = 0, vjust = 0)

# export
ggsave( plot = fd_fig, 
        filename = here::here("5.results_and_figures","figures","functional_diversity.png"),
        width = 16, height = 6, 
        bg = "white",
        dpi=600)
