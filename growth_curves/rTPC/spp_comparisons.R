library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")

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
pairs(d_params.t %>% select(starts_with("ctmax")))
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
########################################################################
########################Completed data set up###########################
########################################################################

# comparisons
# ## max GR ~ MAT, growing season GDD, nongrowing season GDD
# ## Topt ~ MAT, GDD ests
# ## T low ~ Tmin, nongrowing GDD
# ## T high ~ Tmax, growing GDD
# ## breadth50 ~ Tmax-Tmin (intra-annual), stdev of MAT andor GDD (inter-annual variation)
# ## skew ~ Tmax, MAT - Tmax
# ## thermal safety (Tmax-Topt) ~ MAT - Tmax

nrow(d_params.metadata.sp)
colnames(d_params.metadata.sp)

# first compare params ~ spp
# 
# 
# rmax
mod1 = aov(rmax ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#spp          1 0.5786  0.5786   24.49 4.38e-06 ***
#Residuals   76 1.7956  0.0236                     
model.tables(mod1, "means")
#
#        Nd     Nf
#     2.651  2.442
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = rmax
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(paste("Maximum growth rate (mm day"^-1,")"))) +
    annotate(geom = "text", label = "*", x = 1.5, y = 3, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/rmax_spp.pdf")
p1
dev.off()
#
#
#
# topt
mod1 = aov(topt ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value Pr(>F)
#spp          1  0.085  0.0846   0.263   0.61
#Residuals   76 24.445  0.3216                    
model.tables(mod1, "means")
#
#        Nd     Nf
#     21.96  22.04
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = topt
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(paste("Optimum growth temperature (",degree~C,")"))) +
    #annotate(geom = "text", label = "*", x = 1.5, y = 24, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/topt_spp.pdf")
p1
dev.off()

#
#
#
# ctmax.ssh
mod1 = aov(ctmax.ssh ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#spp          1  590.3   590.3   31.36 3.27e-07 ***
#Residuals   76 1430.7    18.8                     

model.tables(mod1, "means")
#
#        Nd     Nf
#     34.9  41.56
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = ctmax.ssh
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(paste("Maximum growth temperature (",degree~C,")"))) +
    annotate(geom = "text", label = "*", x = 1.5, y = 53, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/ctmax.ssh_spp.pdf")
p1
dev.off()
#
#
#
# ctmax.boat
mod1 = aov(ctmax.boat ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#spp          1  109.4  109.43   34.82 9.44e-08 ***
#Residuals   76  238.8    3.14                  

model.tables(mod1, "means")
#
#        Nd     Nf
#     34.34  37.21
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = ctmax.boat
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(atop("Boatman model",paste("Maximum growth temperature (",degree~C,")")))) +
    annotate(geom = "text", label = "*", x = 1.5, y = 41, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/ctmax.boat_spp.pdf")
p1
dev.off()
#
#
#
# ctmax.wt_mean
mod1 = aov(ctmax.wt_mean ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#spp          1  586.7   586.7   31.44 3.17e-07 ***
#Residuals   76 1418.3    18.7                        

model.tables(mod1, "means")
#
#        Nd     Nf
#     34.9  41.54
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = ctmax.wt_mean
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(atop("Weighted mean of best fit models",paste("Maximum growth temperature (",degree~C,")")))) +
    annotate(geom = "text", label = "*", x = 1.5, y = 52, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/ctmax.wt_mean_spp.pdf")
p1
dev.off()

colnames(d_params.metadata.sp)
#
#
#
# t50low
mod1 = aov(t50low ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value Pr(>F)
#spp          1   0.21  0.2141   0.292   0.59
#Residuals   76  55.67  0.7326                           

model.tables(mod1, "means")
#
#        Nd     Nf
#     11.2  11.33
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = t50low
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(atop(paste("Low temperature (",degree~C,")"), "at half maximum growth rate", ))) +
    #annotate(geom = "text", label = "*", x = 1.5, y = 52, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/t50low_spp.pdf")
p1
dev.off()
#
#
#
# t50high
mod1 = aov(t50high ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#spp          1  58.63   58.63   33.14 1.71e-07 ***
#Residuals   76 134.46    1.77                               

model.tables(mod1, "means")
#
#        Nd     Nf
#     27.64  29.74
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = t50high
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(atop(paste("High temperature (",degree~C,")"), "at half maximum growth rate", ))) +
    annotate(geom = "text", label = "*", x = 1.5, y = 32.5, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/t50high_spp.pdf")
p1
dev.off()

colnames(d_params.metadata.sp)
#
#
#
# breadth50
mod1 = aov(breadth50 ~ spp, data = d_params.metadata.sp)
plot(residuals(mod1))
qqnorm(residuals(mod1))
summary(mod1)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#spp          1  51.76   51.76   17.81 6.67e-05 ***
#Residuals   76 220.83    2.91                                    

model.tables(mod1, "means")
#
#        Nd     Nf
#     16.44  18.41
#rep 17.000 61.000
p1 = ggplot(
    d_params.metadata.sp, 
    aes(x = factor(
        spp, 
        levels = c("Nf", "Nd"),
        labels = c("N. faginata", "N. ditissima")
        ), 
        y = breadth50
    )
) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.25) +
    labs(x = "Species", y = expression(paste("Niche breadth (",degree~C,")"))) +
    annotate(geom = "text", label = "*", x = 1.5, y = 24, size = 10) +
    my_gg_theme +
    theme(
        axis.title.x = element_blank()
    )
p1

pdf("figures/rTPC_comparisons/breadth50_spp.pdf")
p1
dev.off()

colnames(d_params.metadata.sp)



