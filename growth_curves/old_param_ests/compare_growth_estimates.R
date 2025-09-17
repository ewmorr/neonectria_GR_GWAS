library(dplyr)
library(ggplot2)
source("R_scripts/ggplot_theme.txt")

#make metadata
source("growth_curves/make_site_metadata.r")
#read gr data
fits.df = read.csv("data/summary_tables/skew_normal_fit_parameters_all_samples.csv")
#assign site
fits.df$Site = sapply(strsplit(fits.df$iso_name, " "), "[[", 1)
#join the gr stats to site metadata (by site)
fits.df.metadata = left_join(
    fits.df,
    site_metadata,
    by = "Site"
)
nrow(fits.df)
nrow(fits.df.metadata)

# join to sample level metadata. Just the species in this case
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[c(3, 5)] = c("iso_name", "spp")
sample_ID_map %>% select(iso_name, spp) %>% unique

fits.df.metadata.sp = left_join(
    fits.df.metadata,
    sample_ID_map %>% select(iso_name, spp) %>% unique,
    by = "iso_name"
)
nrow(fits.df.metadata.sp)
######################
#plots and stats
colnames(fits.df.metadata.sp)


#####################
#maxGR
#####################
#The comnps of maxGR bn spp and with Topt are interesting
# but the comps to heat accumulation are kind of bogus
# what's the hypothesis
#######################

summary(aov(maxGR.exp ~ spp, data = fits.df.metadata.sp))
#spp p = 7.63e-08

p1 = ggplot(fits.df.metadata.sp, aes(spp, maxGR.exp)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    my_gg_theme +
    scale_x_discrete(labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata")) +
    labs(x = "Species", y = expression(paste("Max GR (mm day"^-1, ")" )) ) +
    annotate(geom = "text", label = "*", x = 1.5, y = 3.25, size = 10)
p1

pdf("figures/species_growth_rate.pdf")
p1
dev.off()

#spp*state.name
summary(aov(maxGR.exp ~ spp*state.name, data = fits.df.metadata.sp))
#spp p = 1.16e-07 ***
summary(aov(maxGR.exp ~ state.name, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#p = 0.11

ggplot(fits.df.metadata.sp, aes(state.name, maxGR.exp, color = spp)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    theme_bw()


#Topt
anova(lm(maxGR.exp ~ Topt*spp, data = fits.df.metadata.sp))
#spp p = 10.0002
#Topt p = 5.544e-06
summary(lm(maxGR.exp ~ Topt, data = fits.df.metadata.sp))
#slope p = 0.037, est = 2.36e-05
summary(lm(maxGR.exp ~ Topt, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#p = 1.93e-05, slope = 0.385

p1 = ggplot(fits.df.metadata.sp, aes(Topt, maxGR.exp)) +
    geom_smooth(method = "lm") +
    geom_point(aes(color = spp)) +
    scale_color_manual(values = cbPalette, labels = c("Nd" = "N.ditissima", "Nf" = "N.faginata")) +
    my_gg_theme +
    labs(
        y = expression(paste("Max GR (mm day"^-1, ")" )), 
        x = expression(paste("Growth T"[opt]," (",degree, "C)" ))  
    ) +
    annotate(
        geom = "text", 
        label = "y = 0.01x + 0.21", 
        x = 22.5, y = 0.38, size = 5)+
    annotate(
        geom = "text", 
        label = expression(paste("r"^2, " = 0.1, ", italic(P), " = 0.005")), 
        x = 22.5, y = 0.36, size = 5)
p1

#MAT
anova(lm(maxGR.exp ~ MAT*spp, data = fits.df.metadata.sp))
#spp p = 0.001
summary(lm(maxGR.exp ~ spp/MAT - 1, data = fits.df.metadata.sp))
#no sig slope from 0

ggplot(fits.df.metadata.sp, aes(MAT, maxGR.exp, color = spp)) +
    geom_point() +
    geom_smooth(method = "lm") +
    my_gg_theme +
    scale_color_manual(values = cbPalette)

#GDD-NG
anova(lm(maxGR.exp ~ HDD4.mean_nongrowing*spp, data = fits.df.metadata.sp))
summary(lm(maxGR.exp ~ HDD4.mean_nongrowing*spp, data = fits.df.metadata.sp))
#spp p = 0.0004
summary(lm(maxGR.exp ~ spp/HDD4.mean_nongrowing - 1, data = fits.df.metadata.sp))
#Nf slope marginal sig

ggplot(fits.df.metadata.sp, aes(HDD4.mean_nongrowing, maxGR.exp, color = spp)) +
    geom_point() +
    geom_smooth(method = "lm") +
    my_gg_theme +
    scale_color_manual(values = cbPalette) +
    labs(x = "GDD nongrowing", y = expression(paste("Max GR (mm day"^-1, ")" )) )

#GDD-GS
anova(lm(maxGR.exp ~ HDD4.mean_growing*spp, data = fits.df.metadata.sp))
summary(lm(maxGR.exp ~ HDD4.mean_growing*spp, data = fits.df.metadata.sp))
#spp p = 0.001
summary(lm(maxGR.exp ~ spp/HDD4.mean_growing - 1, data = fits.df.metadata.sp))
#Nf slope marginal sig

ggplot(fits.df.metadata.sp, aes(HDD4.mean_growing, maxGR.exp, color = spp)) +
    geom_point() +
    geom_smooth(method = "lm") +
    my_gg_theme +
    scale_color_manual(values = cbPalette) +
    labs(x = "GDD growing", y = expression(paste("Max GR (mm day"^-1, ")" )) )

##########################
#Topt
###########################

summary(aov(Topt ~ spp, data = fits.df.metadata.sp))
#spp p = 0.006

p1 = ggplot(fits.df.metadata.sp, aes(spp, Topt)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    my_gg_theme +
    scale_x_discrete(labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata")) +
    labs(x = "Species", y = expression(paste("Growth T"[opt]," (",degree, "C)" ))  ) +
    annotate(geom = "text", label = "*", x = 1.5, y = 23.75, size = 10)
p1

pdf("figures/species_Topt.pdf")
p1
dev.off()

#spp*state.name
summary(aov(Topt ~ spp*state.name, data = fits.df.metadata.sp))
#spp p = 0.0045
summary(aov(Topt ~ state.name, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#p = 0.0927

ggplot(fits.df.metadata.sp, aes(state.name, Topt, color = spp)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    theme_bw()

#MAT
mod = lm(Topt ~ MAT*spp, data = fits.df.metadata.sp)
plot(residuals(mod))
qqnorm(residuals(mod))

anova(lm(Topt ~ MAT*spp, data = fits.df.metadata.sp))
#spp p = 0.005
anova(lm(Topt ~ MAT, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
summary(lm(Topt ~ MAT, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#Nf p = 0.048, r-sq = 0.07
anova(lm(Topt ~ MAT, data = fits.df.metadata.sp %>% filter(spp == "Nd")))


p1 = ggplot(fits.df.metadata.sp %>% filter(spp == "Nf"), aes(MAT, Topt)) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    my_gg_theme +
    labs(x = expression(paste("MAT (",degree, "C)" )), 
         y = expression(paste("Growth T"[opt]," (",degree, "C)" )) 
    )  +
    annotate(
        geom = "text", 
        label = "y = 0.11x + 20.0", 
        x = 5, y = 23.1, size = 5)+
    annotate(
        geom = "text", 
        label = expression(paste("r"^2, " = 0.07, ", italic(P), " = 0.048")), 
        x = 5, y = 22.8, size = 5)
p1

pdf("figures/Topt_MAT_Nf.pdf", width = 8, height = 6)
p1
dev.off()


#GDD-NG
summary(lm(Topt ~ HDD4.mean_nongrowing, data = fits.df.metadata.sp))
#ns
summary(lm(Topt ~ HDD4.mean_nongrowing, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#ns
summary(lm(Topt ~ HDD4.mean_nongrowing, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#ns

#GDD-GS
summary(lm(Topt ~ HDD4.mean_growing, data = fits.df.metadata.sp))
#gdd p = 0.09

summary(lm(Topt ~ HDD4.mean_growing, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
mod = lm(Topt ~ HDD4.mean_growing, data = fits.df.metadata.sp %>% filter(spp == "Nf"))
plot(residuals(mod))
qqnorm(residuals(mod))
#p = 0.043, r-sq = 0.068
summary(lm(Topt ~ HDD4.mean_growing, data = fits.df.metadata.sp %>% filter(spp == "Nd")))
#ns


p1 = ggplot(fits.df.metadata.sp %>% filter(spp == "Nf"), aes(HDD4.mean_growing, Topt)) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    my_gg_theme +
    scale_color_manual(values = cbPalette, labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata"), guide = "none") +
    labs(
        x = expression(paste("GDD growing season (",degree, "C)" )) , 
        y = expression(paste("Growth T"[opt], " (",degree, "C)" )) 
    ) +
    annotate(
        geom = "text", 
        label = "y = 0.0006x + 2.0", 
        x = 2000, y = 23.1, size = 5)+
    annotate(
        geom = "text", 
        label = expression(paste("r"^2, " = 0.07, ", italic(P), " = 0.043")), 
        x = 2000, y = 22.8, size = 5)
p1

pdf("figures/Topt_growing_GDD_Nf.pdf", width = 8, height = 6)
p1
dev.off()

###########################
#niche breadth
###########################

summary(aov(niche_breadth ~ spp, data = fits.df.metadata.sp))
# p = 1.27e-06
summary(aov(niche_breadth ~ spp*state.name, data = fits.df.metadata.sp))

p1 = ggplot(fits.df.metadata.sp, aes(spp, niche_breadth)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    my_gg_theme +
    scale_x_discrete(labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata")) +
    labs(x = "Species", y = expression(paste("Niche breadth (",degree, "C)" )) ) +
    annotate(geom = "text", label = "*", x = 1.5, y = 26, size = 10)
p1

pdf("figures/species_niche_breadth.pdf")
p1
dev.off()


#correlation with st.dev. NG
summary(lm(niche_breadth ~ HDD4.sd_nongrowing, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
summary(lm(niche_breadth ~ HDD4.sd_nongrowing, data = fits.df.metadata.sp %>% filter(spp == "Nd")))

#correlation with st.dev. growing
summary(lm(niche_breadth ~ HDD4.sd_growing, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
summary(lm(niche_breadth ~ HDD4.sd_growing, data = fits.df.metadata.sp %>% filter(spp == "Nd")))

#correlation with site temp range
summary(lm(niche_breadth ~ I(tmax-tmin), data = fits.df.metadata.sp %>% filter(spp == "Nf")))
summary(lm(niche_breadth ~ I(tmax-tmin), data = fits.df.metadata.sp %>% filter(spp == "Nd")))


###########################
#site level st dev in Topt
###########################

Topt.sd = fits.df.metadata.sp %>%
    group_by(state.name, spp) %>%
    summarize(Topt.sd = sd(Topt), Topt.mean = mean(Topt))
Topt.sd.metadata = left_join(
    Topt.sd,
    site_metadata[c(-13,-14,-15,-21,-22,-23,-24),],
    by = "state.name"
)

summary(aov(Topt.sd ~ spp, data = Topt.sd.metadata))
# p = 0.02 #this is probably not meaningful given the difs in sampling effort
summary(aov(Topt.mean ~ MAT, data = Topt.sd.metadata %>% filter(spp == "Nf")))

#correlation with st.dev. NG
summary(lm(Topt.sd ~ HDD4.sd_nongrowing, data = Topt.sd.metadata %>% filter(spp == "Nf")))
summary(lm(Topt.sd ~ HDD4.sd_nongrowing, data = Topt.sd.metadata %>% filter(spp == "Nd")))

#correlation with st.dev. growing
summary(lm(Topt.sd ~ HDD4.sd_growing, data = Topt.sd.metadata %>% filter(spp == "Nf")))
summary(lm(Topt.sd ~ HDD4.sd_growing, data = Topt.sd.metadata %>% filter(spp == "Nd")))

#correlation with site temp range
summary(lm(Topt.sd ~ I(tmax-tmin), data = Topt.sd.metadata %>% filter(spp == "Nf")))
summary(lm(Topt.sd ~ I(tmax-tmin), data = Topt.sd.metadata %>% filter(spp == "Nd")))


###########################
#Tlow tmin
###########################

summary(aov(Tlow ~ spp, data = fits.df.metadata.sp))
# ns
summary(aov(Tlow ~ spp*state.name, data = fits.df.metadata.sp))

p1 = ggplot(fits.df.metadata.sp, aes(spp, Tlow)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    my_gg_theme +
    scale_x_discrete(labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata")) +
    labs(x = "Species", y = expression(paste("Growth T"[min]," (",degree, "C)" )) )
p1

pdf("figures/species_Tlow.pdf")
p1
dev.off()


#correlation with tmin
summary(lm(Tlow ~ tmin, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#p = 0.088
summary(lm(Tlow ~ tmin, data = fits.df.metadata.sp %>% filter(spp == "Nd")))
#ns


###########################
#Thigh tmax
###########################

summary(aov(Thigh ~ spp, data = fits.df.metadata.sp))
# 3.36e-08 ***
summary(aov(Thigh ~ spp*state.name, data = fits.df.metadata.sp))

p1 = ggplot(fits.df.metadata.sp, aes(spp, Thigh)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    my_gg_theme +
    scale_x_discrete(labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata")) +
    labs(x = "Species", y = expression(paste("Growth T"[max], " (",degree, "C)" )) ) +
    annotate(geom = "text", label = "*", x = 1.5, y = 35, size = 10)
p1

pdf("figures/species_Thigh.pdf")
p1
dev.off()


#correlation with tmax
summary(lm(Thigh ~ tmax, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#ns
summary(lm(Thigh ~ tmax, data = fits.df.metadata.sp %>% filter(spp == "Nd")))
#ns

###########################
#skewness (k) tmax
###########################

summary(aov(k ~ spp, data = fits.df.metadata.sp))
# 2.8e-08 ***
summary(aov(k ~ spp*state.name, data = fits.df.metadata.sp))
#spp:state.name p = 0.0099

p1 = ggplot(fits.df.metadata.sp, aes(spp, k)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    my_gg_theme +
    scale_x_discrete(labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata")) +
    labs(x = "Species", y = "Skewness (k)" ) +
    annotate(geom = "text", label = "*", x = 1.5, y = 0.55, size = 10)
p1

pdf("figures/species_k.pdf")
p1
dev.off()

p1 = ggplot(fits.df.metadata.sp, aes(state.name, k, color = spp)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    my_gg_theme +
    scale_x_discrete(labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata")) +
    labs(x = "Species", y = expression(paste("Growth T"[max], " (",degree, "C)" )) ) +
    scale_color_brewer(palette = "Paired")
p1

pdf("figures/species_sites_k.pdf", width = 8, height = 6)
p1
dev.off()


#correlation with tmax
summary(lm(k ~ tmax, data = fits.df.metadata.sp %>% filter(spp == "Nf")))
#p = 0.087
summary(lm(k ~ tmax, data = fits.df.metadata.sp %>% filter(spp == "Nd")))
#ns
