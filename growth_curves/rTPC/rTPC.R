library(rTPC)
library(nls.multstart)
library(broom)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)
library(ggrepel)

# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

get_model_names()
?calc_params #explanation of parameters

# preferably include st max, ct min, (britical upper and lower thresholds for growth) topt, and rmax (opt temp and max gr)
? get_thermaltolerance() # Thermal tolerance is calculated as: CTmax - CTmin
? get_breadth() # Thermal performance breadth is calculated as the range of temperatures over which a curve's rate is at least 0.8 of peak. default 0.8 can be changed
? get_skewness() # Skewness is calculated from the values of activation energy (e) and deactivation energy (eh) as: skewness = e - eh

# other potential interesting values are 'e', 'q10'
# 
# ratkowsky_1983    ct min, ct max
# oneill_1972   rmax, ctmax, topt, q10
# modifiedgaussian_2006 rmax,topt
# kamykowski_1985   tmin, tmax (not clear if this one is actually avlbl
# joehnk_2008   rmax,topt
# gaussian_1987 rmax,topt
# briere2_1999    ct min, ct max
# boatman_2017  rmax, tmin, tmax
# 

data("chlorella_tpc")
d <- chlorella_tpc
d
d <- filter(chlorella_tpc, curve_id == 1)
tidyr::nest(d, data = c(temp, rate)) 

###################
#data

#Read the metadata
source("growth_curves/make_site_metadata.r")
#Read the GR data
gr_stat.df = read.table(file = "data/summary_tables/growth_rate_stats.txt", header = T, sep = "\t")
nrow(gr_stat.df)
gr_stat.df$Site = sapply(strsplit(gr_stat.df$iso_name, " "), "[[", 1)
#join the gr stats to metadata
gr_stat.df.metadata = left_join(
    gr_stat.df,
    site_metadata,
    by = "Site"
)

apply(gr_stat.df.metadata, 1, FUN=function(x) all(is.na(x)) ) %>% sum

# join to sample level metadata. Just the species in this case
# sample metadata
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[c(3, 5)] = c("iso_name", "spp")
gr_stat.df.metadata.sp = left_join(
    gr_stat.df.metadata,
    sample_ID_map %>% select(iso_name, spp) %>% unique,
    by = "iso_name"
)
nrow(gr_stat.df.metadata.sp)
# complete.cases is tossing way more so we filter on gr.est specifically
nrow(gr_stat.df.metadata.sp[!is.na(gr_stat.df.metadata.sp$gr.est),])
gr_stat.df.metadata.sp = gr_stat.df.metadata.sp[!is.na(gr_stat.df.metadata.sp$gr.est),]

gr_stat.df.metadata %>% filter(is.na(tmin)) %>% pull(iso_name) %>% unique # need to add climate data for bart and ccm
################################################################################
################################################################################


# test data
d = gr_stat.df.metadata.sp %>% filter(iso_name == "MES1 31.1.1") %>% select(iso_name, gr.est, temp) %>% rename(rate = gr.est)
d


ggplot(d, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')


# fit multiple models
# ratkowsky_1983    ct min, ct max
# oneill_1972   rmax, ctmax, topt, q10
# modifiedgaussian_2006 rmax,topt
# kamykowski_1985   tmin, tmax (not clear if this one is actually avlbl
# joehnk_2008   rmax,topt
# gaussian_1987 rmax,topt
# briere2_1999    ct min, ct max
# boatman_2017  rmax, tmin, tmax

ssfull_start = get_start_vals(d$temp, d$rate, model_name = 'sharpeschoolfull_1981')
sslow_start = get_start_vals(d$temp, d$rate, model_name = 'sharpeschoollow_1981')
get_start_vals(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')
ssfull_start
sslow_start
# the el value is not estimated. 
ssfull_start['el'] = 0
sslow_start['el'] = 0

d_fits <- tidyr::nest(d, data = c(temp, rate)) %>%
  mutate(
      sharpeschoolfull = map(data, ~nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
                        data = .x,
                        iter = c(4,4,4,4,4,4),
                        #start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') - 10,
                        #start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') + 10,
                        start_lower = ssfull_start - 10,
                        start_upper = ssfull_start + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
    sharpeschoollow = map(data, ~nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = sslow_start - 10,
                        start_upper = sslow_start + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
            sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
  )

d_fits.minmax <- tidyr::nest(d, data = c(temp, rate)) %>%
  mutate(
      boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b), #mimax
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 1,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 1,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
      briere2 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),#mimax
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
      kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c), #mimax
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
      ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b), #mimax
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
  )

d_fits <- tidyr::nest(d, data = c(temp, rate)) %>%
  mutate(
      beta = map(data, ~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                        data = .x,
                        iter = c(6,6,6,6,6),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b), #mimax
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 1,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 1,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
      briere2 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),#mimax
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
    joehnk = map(data, ~nls_multstart(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                        data = .x,
                        iter = c(4,4,4,4, 4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') - 1,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') + 1,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
      kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c), #mimax
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
      modifiedgaussian = map(data, ~nls_multstart(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
     oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') ,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') ,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
      ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b), #mimax
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
        quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                        data = .x,
                        iter = c(4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
  )


# look
glimpse(select(d_fits, 1:5))


# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', cols = -iso_name)

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

# get predictions using augment
newdata <- tibble(temp = seq(-5, 40, length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# plot
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = '') +
  geom_hline(aes(yintercept = 0), linetype = 2)


#######################
# model selection
#
#rerun prediction with temp range
# get predictions using augment
# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', -iso_name)

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

# get model fit stats
d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

d_ic

# filter for best model
best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
best_model

# get colour code
col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# plot
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures',
       subtitle= 'The Sharpe-Schoolfield model is the best model') +
  geom_hline(aes(yintercept = 0), linetype = 2) 


calc_params(d_fits$sharpeschoolhigh[[1]])








#rerun prediction with temp range
# get predictions using augment
# stack models
d_stack <- select(d_fits.minmax, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', -iso_name)

# get predictions using augment
newdata <- tibble(temp = seq(-5, 42, length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

# get model fit stats
d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

d_ic

# filter for best model
best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
best_model

# get colour code
col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# plot
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures',
       subtitle= 'The Sharpe-Schoolfield model is the best model') +
  geom_hline(aes(yintercept = 0), linetype = 2) 


calc_params(d_fits.minmax$ratkowsky[[1]])


