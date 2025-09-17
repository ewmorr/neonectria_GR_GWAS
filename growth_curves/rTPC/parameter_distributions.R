library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")
library(diptest)
library(cowplot)


#!!!!!!!!!!!!!!!!!!!!!!
# NOTE that according to ?rTPC::get_ct[min|max] the Sharpe-Schoolfied model 
# min and max estimates are scaled to 5% of the maximum because the model
# never intercepts the x. We should use Boatman for the max. Min is zero
# HOWEVER the SSH model still fits a lot better in some cases especially in 
# cases with heavy skew and is probably still giving better estimates here.
# Maybe this is a case for taking the model avg


###################
#data

#Read the metadata
source("growth_curves/make_site_metadata.r")
#Read the GR data
d_params = read.csv("data/rTPC/all_fits.params.csv")
d_params = d_params %>% filter(model_name %in% c("sharpeschoolhigh", "boatman"))
d_half_max = read.csv("data/rTPC/sharpeschoolhigh_half_max_temps.csv")
colnames(d_half_max)[c(2,4)] = c("t50low", "t50high")
# information criteria with weights
d_ic = read.csv("data/rTPC/d_ic.bestTwo_fits.fit_stats.csv")

iso_names = unique(d_params$iso_name)
d_params.t = data.frame(
    iso_name = vector(mode = "character", length = length(iso_names)),
    rmax = vector(mode = "numeric", length = length(iso_names)),
    topt = vector(mode = "numeric", length = length(iso_names)),
    ctmax.ssh = vector(mode = "numeric", length = length(iso_names)),
    ctmax.boat = vector(mode = "numeric", length = length(iso_names)),
    ctmax.wt_mean = vector(mode = "numeric", length = length(iso_names)),
    breadth = vector(mode = "numeric", length = length(iso_names)),
    skewness = vector(mode = "numeric", length = length(iso_names))
)
for(i in 1:length(iso_names)){
    d_params.t$iso_name[i] = d_params %>% 
        filter(model_name == "sharpeschoolhigh" & iso_name == iso_names[i]) %>% 
        pull(iso_name)
    d_params.t$rmax[i] = d_params %>% 
        filter(model_name == "sharpeschoolhigh" & iso_name == iso_names[i]) %>% 
        pull(rmax)
    d_params.t$topt[i] = d_params %>% 
        filter(model_name == "sharpeschoolhigh" & iso_name == iso_names[i]) %>% 
        pull(topt)
    d_params.t$ctmax.ssh[i] = d_params %>% 
        filter(model_name == "sharpeschoolhigh" & iso_name == iso_names[i]) %>% 
        pull(ctmax)
    d_params.t$ctmax.boat[i] = d_params %>% 
        filter(model_name == "boatman" & iso_name == iso_names[i]) %>% 
        pull(ctmax)
    d_params.t$ctmax.wt_mean[i] = weighted.mean(
        d_params %>% 
            filter(iso_name == iso_names[i]) %>% 
            pull(ctmax), 
        d_ic %>% 
            filter(iso_name == iso_names[i]) %>% 
            pull(weight)
        )
    d_params.t$breadth[i] = d_params %>% 
        filter(model_name == "sharpeschoolhigh" & iso_name == iso_names[i]) %>% 
        pull(breadth)
    d_params.t$skewness[i] = d_params %>% 
        filter(model_name == "sharpeschoolhigh" & iso_name == iso_names[i]) %>% 
        pull(skewness)
}
head(d_params.t)

d_params.t


d_params.t %>% filter(iso_name == "BART C8 2C 287")
d_params.t %>% filter(iso_name == "BART C8 3N 44")
d_params.t %>% filter(iso_name == "CCM ACEPE 10")
d_params.t %>% filter(ctmax.ssh > 50)

d_params = left_join(d_params.t, d_half_max %>% select(-Topt))
head(d_params)
d_params$breadth50 = d_params$t50high - d_params$t50low
head(d_params)

d_params$Site = sapply(strsplit(d_params$iso_name, " "), "[[", 1)
d_params$Site[d_params$Site == "ADS2"] = "ADS1"
d_params$Site[d_params$Site == "MEN2"] = "MEN1"
d_params$Site[d_params$Site == "MES2"] = "MES1"
d_params$Site[d_params$Site == "ADN2"] = "ADN1"
d_params$Site[d_params$Site == "MM"] = "CCM"
#join the gr stats to metadata

d_params.metadata = left_join(
    d_params,
    site_metadata,
    by = "Site"
) 
    

apply(d_params.metadata, 1, FUN=function(x) all(is.na(x)) ) %>% sum

# join to sample level metadata. Just the species in this case
# sample metadata
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[c(3, 5)] = c("iso_name", "spp")
d_params.metadata.sp = left_join(
    d_params.metadata,
    sample_ID_map %>% select(iso_name, spp) %>% unique,
    by = "iso_name"
)
head(d_params.metadata.sp)
sum(complete.cases(d_params.metadata.sp)) == nrow(d_params.metadata.sp)
#d_params.metadata.sp[which(!complete.cases(d_params.metadata.sp)),]

########################################################################
########################Completed data set up###########################
########################################################################

# rmax
# topt
# ctmax.ssh
# ctmax.boat
# breadth
# skewness
# t50low
# t50high
# breadth50

# rmax

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = rmax)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("Maximum growth rate (mm day"^-1,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(rmax))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "rmax")
p2

p_final.rmax <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  


print(p_final.rmax)


shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "rmax"])
# W = 0.96647, p-value = 0.09283
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "rmax"])
# D = 0.030.0392046932, p-value = 0.7267


# topt

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = topt)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("Optimum growth temperature (",degree~C,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(topt))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "Topt")
p2

p_final <- p1 +
  annotation_custom(
    grob = ggplotGrob(p2),
    xmin = 20.6, xmax = 21.25,  # X-range for placement
    ymin = 0.45, ymax = 0.7  # Y-range for placement
  )

print(p_final)

# Combine with cowplot
p_final.topt <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.topt)


shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "topt"])
# W = 0.96217, p-value = 0.05657
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "topt"])
# D = 0.036932, p-value = 0.8154

########################
# ctmax.ssh

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = ctmax.ssh)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("Maximum growth temperature (",degree~C,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(ctmax.ssh))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "ctmax")
p2

# Combine with cowplot
p_final.ctmax <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.ctmax)

shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "ctmax.ssh"])
# W = 0.94079, p-value = 0.005379
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "ctmax.ssh"])
# D = 0.043903, p-value = 0.5251

########################
# ctmax.boat

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = ctmax.boat)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("Boatman maximum growth temperature (",degree~C,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(ctmax.boat))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "ctmax.boat")
p2

# Combine with cowplot
p_final.ctmax.boat <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.ctmax.boat)


shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "ctmax.boat"])
# W = 0.89794, p-value = 9.755e-05
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "ctmax.boat"])
# D = 0.037154, p-value = 0.8081


########################
# breadth

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = breadth)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("Niche breadth at 80% GR (",degree~C,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(breadth))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "breadth80")
p2

# Combine with cowplot
p_final.breadth <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.breadth)


shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "breadth"])
# W = 0.93573, p-value = 0.003192
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "breadth"])
# D = 0.034836, p-value = 0.8844


########################
# skewness

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = skewness)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = "Skewness"
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(skewness))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "skewness")
p2

# Combine with cowplot
p_final.skewness <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.skewness)


shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "skewness"])
# W = 0.64898, p-value = 8.271e-11
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "skewness"])
# D = 0.026634, p-value = 0.9923


########################
# t50low

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = t50low)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("Low temperature at 50% GR (",degree~C,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(t50low))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "t50low")
p2

# Combine with cowplot
p_final.t50low <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.t50low)


shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "t50low"])
# W = 0.97754, p-value = 0.3232
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "t50low"])
# D = 0.037058, p-value = 0.8112


########################
# t50high

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = t50high)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("High temperature at 50% GR (",degree~C,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(t50high))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "t50high")
p2

# Combine with cowplot
p_final.t50high <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.t50high)

shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "t50high"])
# W = 0.89646, p-value = 8.623e-05
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "t50high"])
# D = 0.052812, p-value = 0.2209


########################
# breadth50

p1 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(x = breadth50)
) +
    geom_density(fill = "tomato", alpha = 0.5) +
    labs(
        x = expression(paste("Niche breadth at 50% GR (",degree~C,")"))
    ) +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp %>% filter(spp == "Nf"),
    aes(sample = scale(breadth50))
) +
    stat_qq() +
    stat_qq_line(color = "black", linetype = "dashed") + # theoretical QQ line
    geom_qq(distribution = stats::qnorm) +
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    #facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme +
    labs(y = "breadth50")
p2

# Combine with cowplot
p_final.breadth50 <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2, x = 0.1, y = 0.65, width = 0.25, height = 0.3)  

print(p_final.breadth50)

shapiro.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "breadth50"])
# W = 0.95262, p-value = 0.01923
dip.test(d_params.metadata.sp[d_params.metadata.sp$spp == "Nf", "breadth50"])
# D = 0.035246, p-value = 0.8709


##################################
##################################

# rmax
# topt
# ctmax.ssh
# ctmax.boat
# breadth
# skewness
# t50low
# t50high
# breadth50

pdf("figures/rTPC/parameter_distributions.pdf", width = 10, height = 6)
print(p_final.topt)
print(p_final.ctmax)
print(p_final.ctmax.boat)
print(p_final.breadth)
print(p_final.skewness)
print(p_final.t50low)
print(p_final.t50high)
print(p_final.breadth50)
dev.off()







