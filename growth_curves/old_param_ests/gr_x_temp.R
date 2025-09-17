library(dplyr)
library(ggplot2)
source("growth_curves/make_site_metadata.r")
source("R_scripts/ggplot_theme.txt")

#######################
#Read data

#makes site_metadata
site_metadata
nrow(site_metadata)
site_metadata$Site %>% unique %>% sort

gr_stat.df = read.table(file = "data/summary_tables/growth_rate_stats.txt", header = T, sep = "\t")
nrow(gr_stat.df)
gr_stat.df$Site = sapply(strsplit(gr_stat.df$iso_name, " "), "[[", 1)
gr_stat.df$Site %>% unique %>% sort
gr_stat.df$Site %>% unique %>% length
#join the gr stats to metadata
gr_stat.df.metadata = left_join(
    gr_stat.df,
    site_metadata,
    by = "Site"
)

#check rows on joined and unique sites
nrow(gr_stat.df.metadata)
colnames(gr_stat.df.metadata)
unique(gr_stat.df.metadata$Site)
unique(gr_stat.df.metadata) %>% nrow
#all is well

#check for all na rows
apply(gr_stat.df.metadata, 1, FUN=function(x){ all(is.na(x)) } ) %>% sum

# join to sample level metadata. Just the species in this case
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
sample_ID_map
sample_ID_map %>% nrow
colnames(sample_ID_map)[c(3, 5)] = c("iso_name", "spp")
sample_ID_map %>% select(iso_name, spp) %>% unique

gr_stat.df.metadata.sp = left_join(
    gr_stat.df.metadata,
    sample_ID_map %>% select(iso_name, spp) %>% unique,
    by = "iso_name"
)
nrow(gr_stat.df.metadata.sp)
is.na(gr_stat.df.metadata.sp$spp) %>% sum

gr_stat.df.metadata.sp %>% filter(is.na(spp)) %>% select(iso_name) %>% unique

####################################
#summary plots and stats

#calculate variability between reps of isolates within temps
#start with plotting each rep separately
#then calc s.e. and mean

#color by isolate origin MAT
ggplot(
    gr_stat.df.metadata.sp,
    aes(y = gr.est, x = temp, shape = as.factor(rep), color = MAT )
) +
    geom_point() +
    facet_wrap(~iso_name) +
    my_gg_theme

#sort by spp and site (temp)
gr_stat.df.metadata.sp %>% arrange(rev(spp)) %>% arrange(MAT)

#color by species
p = ggplot(
    gr_stat.df.metadata.sp %>% arrange(rev(spp)) %>% arrange(MAT),
    aes(y = gr.est, x = temp, shape = as.factor(rep), color = spp )
) +
    geom_point() +
    facet_wrap(spp~iso_name) +
    labs(
        y = expression(paste("Growth rate (mm day"^-1, ")")), 
        x = expression(paste("Growth temperature (",degree,C, ")")),
        shape = "Replicate",
        color = "Species"
        ) +
    scale_color_manual(values = cbPalette) +
    theme_bw() 
p

#average reps per temp to calculate optimal temp

pdf("figures/gr_x_spp.pdf", width = 12, height = 12)
p
dev.off()


#add loess fit
p = ggplot(
    gr_stat.df.metadata.sp %>% arrange(rev(spp)) %>% arrange(MAT),
    aes(y = gr.est, x = temp, color = spp )
) +
    geom_point(aes(shape = as.factor(rep))) +
    geom_smooth() +
    facet_wrap(spp~iso_name) +
    labs(
        y = expression(paste("Growth rate (mm day"^-1, ")")), 
        x = expression(paste("Growth temperature (",degree,C, ")")),
        shape = "Replicate",
        color = "Species"
    ) +
    scale_color_manual(values = cbPalette) +
    theme_bw() 
p

#average reps per temp to calculate optimal temp

pdf("figures/gr_x_spp.loess.pdf", width = 12, height = 12)
p
dev.off()


################################
#no. reps per iso x temp

#check for difs by set to make sure no isos shared bn sets
gr_stat.df.metadata.sp %>%
    filter(!is.na(gr.est)) %>%
    group_by(set, iso_name, temp) %>%
    summarize(n_reps = n()) %>% 
    nrow()

gr_stat.df.metadata.sp %>%
    filter(!is.na(gr.est)) %>%
    group_by(iso_name, temp) %>%
    summarize(n_reps = n()) %>% 
    nrow()
# no isos shared

#There are a couple of isos with no valid measurements at 30 C
gr_stat.df.metadata.sp %>%
    group_by(iso_name, temp) %>%
    filter(all(is.na(gr.est))) %>%
    summarize(n_reps = sum(!is.na(gr.est)))

#Count number of not na gr estimates for what needs to be rerun
gr_stat.df.metadata.sp %>%
    group_by(set, iso_name, temp) %>%
    summarize(n_reps = sum(!is.na(gr.est))) %>%
    filter(
        (n_reps < 3 & set == 1) |
            (n_reps < 2 & set == 2)
    ) %>%
    arrange(temp)


#set 1 isos with less than 3 reps at a temp
#or set 2 isos with less than 2 reps at a temp
low_reps = gr_stat.df.metadata.sp %>%
    group_by(set, iso_name, temp) %>%
    summarize(n_reps = sum(!is.na(gr.est))) %>%
    filter(
        (n_reps < 3 & set == 1) |
        (n_reps < 2 & set == 2)
    ) %>%
    arrange(temp)

low_reps

write.table(low_reps, file = "data/summary_tables/low_reps.03072023.csv", sep = ",", col.names = T, row.names = F, quote = F)
