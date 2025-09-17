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

d = gr_stat.df.metadata.sp %>% select(spp, gr.est, temp) %>% rename(rate = gr.est)
head(d)

#make fake zero data to root at zero
#d = rbind(
#    data.frame(
#        spp = c(rep("Nf", 100), rep("Nd", 100)),
#        rate = rep(0,200),
#        temp = rep(0,200)
#    ),
#    d
#)

head(d)

################################################################################
################################################################################

# fit the chosen model formulation in rTPC
d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(
        boatman = {
             map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b), #mimax
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 1,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 1,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        supp_errors = 'Y',
                        convergence_count = FALSE))
             },
      sharpeschoolhigh = {
             map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
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

# create new list column of for high resolution data
d_preds <- mutate(d_fits, 
    new_data = map(
        data, 
        ~tibble(temp = seq(-5, 45, length.out = (length(-5:45)-1)*4+1))
        #~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))
        )
    ) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', -c(spp,new_data)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(spp, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)
head(d_preds)

# plot
p1 = ggplot(d_preds) +
  geom_line(aes(temp, .fitted)) +
  geom_point(aes(temp, rate), d, alpha = 0.25) +
  facet_grid(model_name~spp) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate') +
    theme_bw() 
p1

pdf("figures/rTPC/all_Nd_v_Nf_two_best_models.pdf", width = 8, height = 6)
p1
dev.off()

saveRDS(d_fits, "data/rTPC/Nf_Nd_nls_fits_tibble.rds")
saveRDS(data.frame(d_preds), "data/rTPC/Nf_Nd_nls_preds_df.rds")

d_params <- pivot_longer(d_fits, names_to = 'model_name', values_to = 'fit', -c(spp,data)) %>%
  mutate(params = map(fit, calc_params)) %>%
  select(spp, model_name, params) %>%
  unnest(params)
d_params

saveRDS(data.frame(d_params %>% select(-q10, thermal_safety_margin)), "data/rTPC/Nf_Nd_nls_params_df.rds")

    