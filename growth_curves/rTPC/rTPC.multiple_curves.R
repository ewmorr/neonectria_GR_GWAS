library(rTPC)
library(nls.multstart)
library(broom)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)
library(ggrepel)
library(magrittr)

# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

# when scaling up our code to fit hundreds of models, its nice to have a progress bar
# edit nls_multstart to allow for a progress bar
nls_multstart_progress <- function(formula, data = parent.frame(), iter, start_lower, 
                                   start_upper, supp_errors = c("Y", "N"), convergence_count = 100, 
                                   control, modelweights, ...){
  if(!is.null(pb)){
    pb$tick()
  }
  nls_multstart(formula = formula, data = data, iter = iter, start_lower = start_lower, 
                start_upper = start_upper, supp_errors = supp_errors, convergence_count = convergence_count, 
                control = control, modelweights = modelweights, ...)
}

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

gr_stat.df.metadata %>% filter(is.na(tmin)) %>% pull(iso_name) %>% unique 

d = gr_stat.df.metadata.sp %>% select(iso_name, gr.est, temp) %>% rename(rate = gr.est)
################################################################################
################################################################################

# start progress bar and estimate time it will take
number_of_models <- 10
number_of_curves <- length(unique(d$i))

# setup progress bar
pb <- progress::progress_bar$new(total = number_of_curves*number_of_models,
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")

# fit the chosen model formulation in rTPC
d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(
      beta = {
          map(data, ~nls_multstart_progress(rate~beta_2012(temp = temp, a, b, c, d, e),
                        data = .x,
                        iter = c(6,6,6,6,6),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
          },
         boatman = {
             map(data, ~nls_multstart_progress(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b), #mimax
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 1,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 1,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
             },
      briere2 = {
          map(data, ~nls_multstart_progress(rate~briere2_1999(temp = temp, tmin, tmax, a,b),#mimax
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
          },
    joehnk = {
        map(data, ~nls_multstart_progress(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                        data = .x,
                        iter = c(4,4,4,4, 4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') - 1,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') + 1,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
        },
      kamykowski = {
          map(data, ~nls_multstart_progress(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c), #mimax
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
          },
      modifiedgaussian = {
          map(data, ~nls_multstart_progress(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
          },
     oneill = {
         map(data, ~nls_multstart_progress(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') ,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') ,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
         },
      ratkowsky = {
          map(data, ~nls_multstart_progress(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b), #mimax
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
          },
        quadratic = {
            map(data, ~nls_multstart_progress(rate~quadratic_2008(temp = temp, a, b, c),
                        data = .x,
                        iter = c(4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
            },
         sharpeschoolhigh = {
             map(data, ~nls_multstart_progress(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
         }
      )
#saveRDS(d_fits, "data/rTPC/all_fits.rds")
d_fits = readRDS("data/rTPC/all_fits.rds")

# also fit the shapreschoolfull
# need to manually set start vals
ssfull_start = get_start_vals(d$temp, d$rate, model_name = 'sharpeschoolfull_1981')
ssfull_start['el'] = 0

d_fits.ssfull <- nest(d, data = c(temp, rate)) %>%
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
                        convergence_count = FALSE))
  )
d_fits$sharpeschoolfull = d_fits.ssfull$sharpeschoolfull

#saveRDS(d_fits, "data/rTPC/all_fits.rds")
d_fits = readRDS("data/rTPC/all_fits.rds")


# create new list column of for high resolution data
d_preds <- mutate(d_fits, 
    new_data = map(
        data, 
        ~tibble(temp = seq(-5, 40, length.out = 100))
        #~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))
        )
    ) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', -c(iso_name,new_data)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(iso_name, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot
p1 = ggplot(d_preds) +
  geom_line(aes(temp, .fitted)) +
  geom_point(aes(temp, rate), d) +
  facet_grid(model_name~iso_name, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'All fitted thermal performance curves') +
    theme_bw() +
    ylim(-0.5,3.5)
p1

isonames = d_preds$iso_name %>% unique
i=1
ggplot() +
    geom_line(
        data = d_preds %>% filter(iso_name == isonames[i]),
        aes(temp, .fitted), 
    ) +
    geom_point(
        data = d %>% filter(iso_name == isonames[i]),
        aes(temp, rate)
    ) +
    facet_wrap(~model_name, ncol = 6) +
    theme_bw() +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = isonames[i]) +
    theme_bw() +
    ylim(-0.5,3.5)



pdf("figures/rTPC/all_rTPCs.pdf", width = 12, height = 8)
for(i in 1:length(isonames)){
print(
    ggplot() +
    geom_line(
        data = d_preds %>% filter(iso_name == isonames[i]),
        aes(temp, .fitted), 
    ) +
    geom_point(
        data = d %>% filter(iso_name == isonames[i]),
        aes(temp, rate)
    ) +
    facet_wrap(~model_name, ncol = 5) +
    theme_bw() +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = isonames[i]) +
    theme_bw() +
    ylim(-0.5,3.5)
)
}
dev.off()

# ratkowsky fails catastrophically for some isos
# # the max temp curve dips then shoots back up (i.e., not unimodal)

# stack models and calculate extra params
d_params <- pivot_longer(d_fits, names_to = 'model_name', values_to = 'fit', -c(iso_name,data)) %>%
  mutate(params = map(fit, calc_params)) %>%
  select(iso_name, model_name, params) %>%
  unnest(params)
write.csv(d_params, "data/rTPC/all_fits.params.csv", row.names = F, quote = F)
##################################################################
# lets also get the temps at 80% of opt (80% is used for breadth)
# Only running sharpeschoolhigh since that is what we will use
# 
# niche breadth is the T range at half max GR
# we will calculate yest before the function
# temp is the input to the predictions function (dense temp series)
# maxGR estimated by max(yest)
# x is the cutoff. input 0.8 for GR at 80% of max
f.breadth_at_x_max = function(temps, yest, x) {
    
    #xMax GR for estimation
    xMax = max(yest)*(1-x)
    
    #break the temp vector and gr vec at the highpoint in GR
    temps.low = temps[ 1:which.max(yest) ] #lower half
    temps.high = temps[ which.max(yest):length(temps) ] #upper half
    yest.low = yest[ 1:which.max(yest) ] #lower half
    yest.high = yest[ which.max(yest):length(temps) ] #upper half
    
    #get the temp at minimum difs between yest and xMax or max based on absolute values
    Topt = temps[ which.max(yest) ]
    Tlow = temps.low[ which.min(abs(xMax - yest.low)) ] 
    Thigh = temps.high[ which.min(abs(xMax - yest.high)) ] 
    
    return(c(Tlow, Topt, Thigh))
}

# create new list column of for high resolution data
d_preds <- mutate(d_fits, 
    new_data = map(
        data, 
        ~tibble(temp = seq(0, 40, length.out = (length(0:40)-1)*100+1))
        )
    ) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', -c(iso_name,new_data)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(iso_name, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)
head(d_preds)

d_preds %>% 
    filter(iso_name == "MES1 31.1.1" & model_name == "sharpeschoolhigh") %$%
    f.breadth_at_x_max(temps = temp, yest = .fitted, x = 0.8)

d_preds %>% 
    filter(model_name == "sharpeschoolhigh") -> d_pred.ssh
iso_names = unique(d_preds$iso_name)
ssh_50T = data.frame(
    iso_name = vector(mode = "character", length = length(iso_names)),
    Tlow = vector(mode = "numeric", length = length(iso_names)),
    Topt = vector(mode = "numeric", length = length(iso_names)),
    Thigh = vector(mode = "numeric", length = length(iso_names))
)

pb <- progress::progress_bar$new(
    total = length(iso_names),
    clear = FALSE,
    format ="[:bar] :percent :elapsedfull"
)

for(i in 1:length(iso_names)){
    d_pred.ssh.t = d_pred.ssh %>% filter(iso_name == iso_names[i])
    ssh_50T$iso_name[i] = iso_names[i]
    Tests = f.breadth_at_x_max(temps = d_pred.ssh.t$temp, yest = d_pred.ssh.t$.fitted, x = 0.5)
    ssh_50T$Tlow[i] = Tests[1]
    ssh_50T$Topt[i] = Tests[2]
    ssh_50T$Thigh[i] = Tests[3]
    pb$tick()
}
head(ssh_50T)

write.csv(ssh_50T, "data/rTPC/sharpeschoolhigh_half_max_temps.csv", quote = F, row.names = F)

minMaxModels = c(
    "boatman",
    "briere2",
    "kamykowski",
    "ratkowsky"
)

# get model fit stats
d_ic <- pivot_longer(d_fits, names_to = 'model_name', values_to = 'fit', -c(iso_name,data)) %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(iso_name, model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic <- d_ic %>%
  # filter(d_ic, aic - min(aic) <= 2) %>%
  mutate(., weight = MuMIn::Weights(AICc))
write.csv(d_ic, "data/rTPC/all_fits.fit_stats.csv", quote = F, row.names = F)

# get model fit stats for JUST SSH and Boat to calculate weight s for just these
d_ic.bestTwo <- pivot_longer(d_fits, names_to = 'model_name', values_to = 'fit', c(boatman,sharpeschoolhigh)) %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(iso_name, model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic.bestTwo <- d_ic.bestTwo %>%
  # filter(d_ic, aic - min(aic) <= 2) %>%
  mutate(., weight = MuMIn::Weights(AICc))
write.csv(d_ic.bestTwo, "data/rTPC/d_ic.bestTwo_fits.fit_stats.csv", quote = F, row.names = F)



d_ic$minmax = ifelse(d_ic$model_name %in% minMaxModels, "Tmin_Tmax", "other")

p2 = d_ic %>%
    pivot_longer(names_to = "fit_stat", values_to = "values", cols = -c(iso_name, model_name, minmax, df.residual)) %>%
    ggplot(.,
        aes(x = model_name, y = values, color = minmax)
    ) +
    facet_wrap(~fit_stat, scales = "free_y", strip.position = "left") +
    geom_boxplot() +
    scale_color_brewer(palette = "Dark2") +
    geom_point(position = position_jitter(width = 0.25), alpha = 0.25, size = 1) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 10)
    )
p2
pdf("figures/rTPC/model_performance.pdf", width = 10, height = 6)
p2
dev.off()

# an example of an isolate where ratkowski fails
d_params %>% filter(iso_name == "ADN1 19.1.1")

pdf("figures/rTPC/best_models_compare_Tmax.pdf", width = 6, height = 4)
d_params %>%
    filter(model_name %in% c("boatman", "sharpeschoolhigh")) %>%
    select(model_name, iso_name, ctmax) %>%
    pivot_wider(names_from = "model_name", values_from = "ctmax") %>%
    ggplot(., aes(x = boatman, y = sharpeschoolhigh)) +
    geom_point() +
    labs(x = "boatman Tmax", y = "shapreschoolhigh Tmax") +
    theme_bw() +
    annotate(
        geom = "text", 
        label = expression(paste("r = 0.96, df = 76, ",italic(P)," < 0.0001")),
        x = 34, y = 51
    )
dev.off()

d_params %>%
    filter(model_name %in% c("boatman", "sharpeschoolhigh")) %>%
    select(model_name, iso_name, ctmax) %>%
    pivot_wider(names_from = "model_name", values_from = "ctmax") -> d_params.best2

    cor(x = d_params.best2$boatman, y = d_params.best2$sharpeschoolhigh)
    cor.test(x = d_params.best2$boatman, y = d_params.best2$sharpeschoolhigh)    

    mod1 = lm(sharpeschoolhigh ~ boatman, data = d_params.best2)
    summary(mod1)

    
    
    
# new preds for plotting schoolhouse up to 50C
    
# create new list column of for high resolution data
d_preds <- mutate(d_fits, 
    new_data = map(
        data, 
        ~tibble(temp = seq(-5, 50, length.out = 110))
        #~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))
        )
    ) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', -c(iso_name,new_data)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(iso_name, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot
p1 = ggplot(d_preds %>% filter(model_name == "sharpeschoolhigh")) +
  geom_line(aes(temp, .fitted)) +
  geom_point(aes(temp, rate), d) +
  facet_wrap(~iso_name) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'All Schoolhouse-high thermal performance curves') +
    theme_bw() +
    geom_hline(yintercept = 0, color = "blue", linetype = 2)
p1

pdf("figures/rTPC/all_schoolhousehigh.pdf", width = 12, height = 10)
p1
dev.off()

# plot
p2 = ggplot(d_preds %>% filter(model_name == "boatman")) +
  geom_line(aes(temp, .fitted)) +
  geom_point(aes(temp, rate), d) +
  facet_wrap(~iso_name) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'All Boatman thermal performance curves') +
    theme_bw() +
    geom_hline(yintercept = 0, color = "blue", linetype = 2)
p2

pdf("figures/rTPC/all_boatman.pdf", width = 12, height = 10)
p2
dev.off()

#############################
# The take home is: sharpeschoolhigh is the best fit (by quite a bit)
# BUT the low temp estimates appear to be garbage
# That is partly a product of "not estimating rate in the area of Arrhenius non-linearity"
# as described in the section of the Schoolhouse paper describing the 4 parameter models
# (note that we also looked at schoolhouse low and the fit on the high side of the curve is garbage)
# BUT if we assume a linear fit we should have 0 in most cases.
# Maybe we take both the Tmax ests from schoolhouse high AND the est from Boatman,
# and use schoolhouse for all of the other parameter ests