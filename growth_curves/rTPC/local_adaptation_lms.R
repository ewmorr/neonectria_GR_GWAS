library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")
library(nlme)

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
#d_params.metadata.sp[which(!complete.cases(d_params.metadata.sp)),]

########################################################################
########################Completed data set up###########################
########################################################################

# comparisons
# ## max GR ~ MAT, growing season GDD, nongrowing season GDD ### This doesn't really make sense (Topt does) but we run anyways
# ## Topt ~ MAT, GDD ests
# ## T low ~ Tmin, nongrowing GDD
# ## T high ~ Tmax, growing GDD
# ## ctmax ~ Tmax
# ## breadth50 ~ Tmax-Tmin (intra-annual), stdev of MAT andor GDD (inter-annual variation)
# ## skew ~ Tmax, MAT - Tmax
# ## thermal safety (Tmax-Topt) ~ MAT - Tmax

nrow(d_params.metadata.sp)
colnames(d_params.metadata.sp)

# rmax

mod1 = lm(rmax ~ MAT*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

mod1 = lm(rmax ~ HDD4.mean_growing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

mod1 = lm(rmax ~ HDD4.mean_nongrowing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
# none sig


#######################################################
# topt

# ~ MAT
mod1 = lm(topt ~ MAT*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Call:
#lm(formula = topt ~ MAT * spp, data = d_params.metadata.sp)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-1.58318 -0.31036  0.07745  0.38612  0.95609 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 22.63365    0.65890  34.350   <2e-16 ***
#MAT         -0.09606    0.09271  -1.036   0.3035    
#sppNf       -1.52788    0.71906  -2.125   0.0369 *  
#MAT:sppNf    0.22130    0.09995   2.214   0.0299 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.5321 on 74 degrees of freedom
#Multiple R-squared:  0.1457,	Adjusted R-squared:  0.1111 
#F-statistic: 4.207 on 3 and 74 DF,  p-value: 0.008363

anova(mod1)
#Analysis of Variance Table
#
#Response: topt
#          Df  Sum Sq Mean Sq F value   Pr(>F)   
#MAT        1  2.1738 2.17376  7.6762 0.007071 **
#spp        1  0.0122 0.01222  0.0431 0.836021   
#MAT:spp    1  1.3883 1.38826  4.9024 0.029900 * 
#Residuals 74 20.9554 0.28318                    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#

mod2 = lme(fixed = topt ~ MAT*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)

# this gives the same sig levels (not the same P of course)
#Linear mixed-effects model fit by REML
#  Data: d_params.metadata.sp 
#       AIC      BIC    logLik
#  144.1005 157.9249 -66.05024
#
#Random effects:
# Formula: ~1 | state.name
#        (Intercept)  Residual
#StdDev:   0.1273674 0.5188116
#
#Fixed effects:  topt ~ MAT * spp 
#                Value Std.Error DF  t-value p-value
#(Intercept) 22.655448 0.6760695 57 33.51053  0.0000
#MAT         -0.100367 0.0952905 17 -1.05327  0.3070
#sppNf       -1.540594 0.7247515 57 -2.12569  0.0379
#MAT:sppNf    0.224020 0.1012487 57  2.21257  0.0309
# Correlation: 
#          (Intr) MAT    sppNf 
#MAT       -0.979              
#sppNf     -0.904  0.888       
#MAT:sppNf  0.896 -0.918 -0.977
#
#Standardized Within-Group Residuals:
#       Min         Q1        Med         Q3        Max 
#-2.9238011 -0.5162359  0.1705788  0.7510737  1.7424936 
#
#Number of Observations: 78
#Number of Groups: 19 

anova(mod2)
#
#            numDF denDF   F-value p-value
#(Intercept)     1    57 110409.74  <.0001
#MAT             1    17      6.25  0.0229
#spp             1    57      0.03  0.8647
#MAT:spp         1    57      4.90  0.0309


p1 = ggplot(d_params.metadata.sp, 
    aes(
        x = MAT, 
        y = topt, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,1), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean annual temperature (",degree~C,")")), 
        y = expression(paste("Optimum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p1


# ~ HDD4.mean_growing
# 
mod1 = lm(topt ~ HDD4.mean_growing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#Call:
#lm(formula = topt ~ HDD4.mean_growing * spp, data = d_params.metadata.sp)

#Residuals:
#    Min      1Q  Median      3Q     Max 
#-1.6231 -0.2980  0.1194  0.3787  0.9968 
#
#Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             23.3895592  1.2093953  19.340   <2e-16 ***
#HDD4.mean_growing       -0.0006355  0.0005361  -1.185   0.2396    
#sppNf                   -2.7894277  1.2950247  -2.154   0.0345 *  
#HDD4.mean_growing:sppNf  0.0012634  0.0005719   2.209   0.0303 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.5352 on 74 degrees of freedom
#Multiple R-squared:  0.1359,	Adjusted R-squared:  0.1009 
#F-statistic:  3.88 on 3 and 74 DF,  p-value: 0.01239

anova(mod1)
#Response: topt
#                      Df  Sum Sq Mean Sq F value  Pr(>F)  
#HDD4.mean_growing      1  1.8988 1.89879  6.6291 0.01203 *
#spp                    1  0.0371 0.03708  0.1294 0.72003  
#HDD4.mean_growing:spp  1  1.3979 1.39788  4.8803 0.03026 *
#Residuals             74 21.1959 0.28643                  


mod2 = lme(fixed = topt ~ HDD4.mean_growing*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#again same results as raw anova


p2 = ggplot(d_params.metadata.sp, 
    aes(
        x = HDD4.mean_growing, 
        y = topt, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,1), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Growing season GDD"[4]," (",degree~C,")")), 
        y = expression(paste("Optimum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p2

# ~ HDD4.mean_nongrowing
# 

mod1 = lm(topt ~ HDD4.mean_nongrowing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#
#Call:
#lm(formula = topt ~ HDD4.mean_nongrowing * spp, data = d_params.metadata.sp)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-1.53079 -0.30534  0.02779  0.44479  1.01636 
#
#Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                22.687730   2.246515  10.099 1.44e-15 ***
#HDD4.mean_nongrowing       -0.005173   0.016033  -0.323    0.748    
#sppNf                      -2.384568   2.335999  -1.021    0.311    
#HDD4.mean_nongrowing:sppNf  0.016869   0.016593   1.017    0.313    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.5468 on 74 degrees of freedom
#Multiple R-squared:  0.09805,	Adjusted R-squared:  0.06149 
#F-statistic: 2.682 on 3 and 74 DF,  p-value: 0.05293
anova(mod1)
#Response: topt
#                         Df  Sum Sq Mean Sq F value  Pr(>F)  
#HDD4.mean_nongrowing      1  2.0350 2.03498  6.8065 0.01098 *
#spp                       1  0.0028 0.00281  0.0094 0.92300  
#HDD4.mean_nongrowing:spp  1  0.3674 0.36745  1.2290 0.27119  
#Residuals                74 22.1244 0.29898                  
#
#

mod2 = lme(fixed = topt ~ HDD4.mean_nongrowing*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
# same


p3 = ggplot(d_params.metadata.sp, 
    aes(
        x = HDD4.mean_nongrowing, 
        y = topt, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Nongrowing season GDD"[4]," (",degree~C,")")), 
        y = expression(paste("Optimum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p3

p4 = ggplot(d_params.metadata.sp,
    aes(x = topt)
) +
    geom_histogram(binwidth = 0.2) +
    labs(
        x = expression(paste("Optimum growth temperature (",degree~C,")"))
    ) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p4

p5 = ggplot(d_params.metadata.sp, 
    aes(
        x = MAT, 
        y = topt, 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_smooth(method = "lm", alpha = 0.25, color = "black") +
    geom_point(aes(color = state.name)) +
    scale_linetype_manual(values = c(2,1), guide = "none") +
    scale_color_manual(values = c25) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
       x = expression(paste("Mean annual temperature (",degree~C,")")), 
        y = expression(paste("Optimum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p5

p6 = ggplot(d_params.metadata.sp, 
    aes(
        x = HDD4.mean_growing, 
        y = topt, 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_smooth(method = "lm", alpha = 0.25, color = "black") +
    geom_point(aes(color = state.name)) +
    scale_linetype_manual(values = c(2,1), guide = "none") +
    scale_color_manual(values = c25) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
       x = expression(paste("Growing season GDD"[4]," (",degree~C,")")), 
        y = expression(paste("Optimum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p6

pdf("figures/rTPC_local_adap_comps/topt.pdf", width = 11, height = 5)
p1
p2
p3
p4
p5
p6
dev.off()
    
###########################################################################
# ## t50low ~ Tmin, nongrowing GDD
# 
# 

mod1 = lm(t50low ~ tmin*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-2.83613 -0.51471  0.05678  0.49702  1.84749 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 10.40660    0.91972  11.315   <2e-16 ***
#tmin        -0.06674    0.07513  -0.888    0.377    
#sppNf        1.46622    0.97944   1.497    0.139    
#tmin:sppNf   0.11890    0.08114   1.465    0.147    
#
#Residual standard error: 0.8465 on 74 degrees of freedom
#Multiple R-squared:  0.05113,	Adjusted R-squared:  0.01266 
#F-statistic: 1.329 on 3 and 74 DF,  p-value: 0.2714
#

mod2 = lme(fixed = t50low ~ tmin*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same


mod1 = lm(t50low ~ HDD4.mean_nongrowing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Residuals:
#    Min      1Q  Median      3Q     Max 
#-2.7398 -0.5314  0.1430  0.4634  1.8568 
#
#Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                17.30188    3.42104   5.057 2.99e-06 ***
#HDD4.mean_nongrowing       -0.04360    0.02441  -1.786   0.0782 .  
#sppNf                      -7.63313    3.55730  -2.146   0.0352 *  
#HDD4.mean_nongrowing:sppNf  0.05476    0.02527   2.167   0.0334 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.8336 on 74 degrees of freedom
#Multiple R-squared:  0.07999,	Adjusted R-squared:  0.0427 
#F-statistic: 2.145 on 3 and 74 DF,  p-value: 0.1018
anova(mod1)
#Response: t50low
#                         Df Sum Sq Mean Sq F value  Pr(>F)  
#HDD4.mean_nongrowing      1  1.163  1.1628  1.6735 0.19982  
#spp                       1  0.044  0.0444  0.0639 0.80117  
#HDD4.mean_nongrowing:spp  1  3.264  3.2636  4.6969 0.03343 *
#Residuals                74 51.418  0.6948                  

mod2 = lme(fixed = t50low ~ HDD4.mean_nongrowing*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#marginal sig on the interaction



p1 = ggplot(d_params.metadata.sp, 
    aes(
        x = tmin, 
        y = t50low, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean minimum temperature of coldest month (",degree~C,")")), 
        y = expression(atop(paste("Low temperature (",degree~C,")"), "at half maximum growth rate", ))
    ) +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp, 
    aes(
        x = HDD4.mean_nongrowing, 
        y = t50low, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,1), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Nongrowing season GDD"[4]," (",degree~C,")")), 
        y = expression(atop(paste("Low temperature (",degree~C,")"), "at half maximum growth rate", ))
    ) +
    my_gg_theme
p2

p4 = ggplot(d_params.metadata.sp, 
    aes(
        x = HDD4.mean_nongrowing, 
        y = t50low, 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_smooth(method = "lm", alpha = 0.25, color = "black") +
    geom_point(aes(color = state.name)) +
    scale_linetype_manual(values = c(2,1), guide = "none") +
    scale_color_manual(values = c25) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
       x = expression(paste("Nongrowing season GDD"[4]," (",degree~C,")")), 
        y = expression(atop(paste("Low temperature (",degree~C,")"), "at half maximum growth rate"))
    ) +
    my_gg_theme
p4

p3 = ggplot(d_params.metadata.sp,
    aes(x = t50low)
) +
    geom_histogram(binwidth = 0.5) +
    labs(
        x = expression(atop(paste("Low temperature (",degree~C,")"), "at half maximum growth rate", ))
    ) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p3

pdf("figures/rTPC_local_adap_comps/t50low.pdf", width = 11, height = 5)
p1
p2
p4
p3
dev.off()


###########################################################################
# ## t50high ~ Tmax, growing GDD
# 
# 

mod1 = lm(t50high ~ tmax*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-4.6054 -0.6695  0.2535  0.9109  2.0688 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept) 25.30737    7.74426   3.268  0.00164 **
#tmax         0.08859    0.29343   0.302  0.76355   
#sppNf        2.97338    8.27586   0.359  0.72041   
#tmax:sppNf  -0.03252    0.31396  -0.104  0.91779   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.345 on 74 degrees of freedom
#Multiple R-squared:  0.3069,	Adjusted R-squared:  0.2788 
#F-statistic: 10.92 on 3 and 74 DF,  p-value: 5.092e-06

mod2 = lme(fixed = t50high ~ tmax*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same


mod1 = lm(t50high ~ HDD4.mean_growing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#Call:
#lm(formula = t50high ~ HDD4.mean_growing * spp, data = d_params.metadata.sp)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-4.8024 -0.6443  0.2303  0.9765  2.3301 
#
#Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             26.4179091  3.0679351   8.611 9.04e-13 ***
#HDD4.mean_growing        0.0005470  0.0013615   0.402    0.689    
#sppNf                    3.9352930  3.2808691   1.199    0.234    
#HDD4.mean_growing:sppNf -0.0008121  0.0014505  -0.560    0.577    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.344 on 74 degrees of freedom
#Multiple R-squared:  0.3078,	Adjusted R-squared:  0.2797 
#F-statistic: 10.97 on 3 and 74 DF,  p-value: 4.853e-06

mod2 = lme(fixed = t50high ~ HDD4.mean_growing*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same


p1 = ggplot(d_params.metadata.sp, 
    aes(
        x = tmax, 
        y = t50high, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean maximum temperature of warmest month (",degree~C,")")), 
        y = expression(atop(paste("High temperature (",degree~C,")"), "at half maximum growth rate", ))
    ) +
    my_gg_theme
p1

p4 = ggplot(d_params.metadata.sp, 
    aes(
        x = tmax, 
        y = t50high, 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_smooth(method = "lm", alpha = 0.25, color = "black") +
    geom_point(aes(color = state.name)) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_manual(values = c25) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean maximum temperature of warmest month (",degree~C,")")), 
        y = expression(atop(paste("High temperature (",degree~C,")"), "at half maximum growth rate", ))
    ) +
    my_gg_theme
p4

p2 = ggplot(d_params.metadata.sp, 
    aes(
        x = HDD4.mean_growing, 
        y = t50high, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Growing season GDD"[4]," (",degree~C,")")), 
        y = expression(atop(paste("High temperature (",degree~C,")"), "at half maximum growth rate", ))
    ) +
    my_gg_theme
p2

p3 = ggplot(d_params.metadata.sp,
    aes(x = t50high)
) +
    geom_histogram(binwidth = 0.5) +
    labs(
        x = expression(atop(paste("High temperature (",degree~C,")"), "at half maximum growth rate", ))
    ) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p3

pdf("figures/rTPC_local_adap_comps/t50high.pdf", width = 11, height = 5)
p1
p4
p2
p3
dev.off()

###########################################################################
# ## ctmax ~ Tmax

mod1 = lm(ctmax.ssh ~ tmax*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Call:
#lm(formula = ctmax.ssh ~ tmax * spp, data = d_params.metadata.sp)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-12.0573  -2.0432   0.7972   2.7862  10.7229 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  21.6567    25.2732   0.857    0.394
#tmax          0.5023     0.9576   0.525    0.601
#sppNf        20.1939    27.0081   0.748    0.457
#tmax:sppNf   -0.5132     1.0246  -0.501    0.618
#
#Residual standard error: 4.389 on 74 degrees of freedom
#Multiple R-squared:  0.2947,	Adjusted R-squared:  0.2661 
#F-statistic: 10.31 on 3 and 74 DF,  p-value: 9.517e-06

mod2 = lme(fixed = ctmax.ssh ~ tmax*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same


mod1 = lm(ctmax.boat ~ tmax*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Call:
#lm(formula = ctmax.boat ~ tmax * spp, data = d_params.metadata.sp)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-5.3661 -0.6966  0.2378  1.0700  2.9226 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  30.5093    10.3352   2.952  0.00423 **
#tmax          0.1454     0.3916   0.371  0.71142   
#sppNf         6.4963    11.0446   0.588  0.55820   
#tmax:sppNf   -0.1375     0.4190  -0.328  0.74376   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.795 on 74 degrees of freedom
#Multiple R-squared:  0.3155,	Adjusted R-squared:  0.2878 
#F-statistic: 11.37 on 3 and 74 DF,  p-value: 3.24e-06

mod2 = lme(fixed = ctmax.boat ~ tmax*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same


mod1 = lm(ctmax.ssh ~ HDD4.mean_growing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)   
#(Intercept)             29.104741   9.768265   2.980   0.0039 **
#HDD4.mean_growing        0.002585   0.004330   0.597   0.5524   
#sppNf                   17.960949  10.459892   1.717   0.0901 . 
#HDD4.mean_growing:sppNf -0.004977   0.004619  -1.077   0.2848   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.323 on 74 degrees of freedom
#Multiple R-squared:  0.3158,	Adjusted R-squared:  0.2881 
#F-statistic: 11.39 on 3 and 74 DF,  p-value: 3.192e-06

mod2 = lme(fixed = ctmax.ssh ~ HDD4.mean_growing*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same

mod1 = lm(ctmax.boat ~ HDD4.mean_growing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             33.4127107  4.0475123   8.255 4.26e-12 ***
#HDD4.mean_growing        0.0004157  0.0017963   0.231    0.818    
#sppNf                    5.8991429  4.3284351   1.363    0.177    
#HDD4.mean_growing:sppNf -0.0013284  0.0019136  -0.694    0.490    

mod1 = lm(ctmax.ssh ~ bio9*spp, data = d_params.metadata.sp) # temp of driest
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  35.6981     2.5247  14.139   <2e-16 ***
#bio9          0.1595     0.4594   0.347   0.7295    
#sppNf         5.4427     2.6055   2.089   0.0402 *  
#bio9:sppNf   -0.3024     0.4722  -0.640   0.5239    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.343 on 74 degrees of freedom
#Multiple R-squared:  0.3094,	Adjusted R-squared:  0.2814 
#F-statistic: 11.05 on 3 and 74 DF,  p-value: 4.475e-06

mod1 = lm(ctmax.ssh ~ bio8*spp, data = d_params.metadata.sp) # temp of wettest
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Call:
#lm(formula = ctmax.ssh ~ bio8 * spp, data = d_params.metadata.sp)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-12.2611  -2.0225   0.8571   2.7323  10.8308 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  35.9164     1.7358  20.692  < 2e-16 ***
#bio8         -0.1334     0.1805  -0.739  0.46216    
#sppNf         5.8949     2.1023   2.804  0.00644 ** 
#bio8:sppNf    0.1122     0.2017   0.557  0.57955    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.379 on 74 degrees of freedom
#Multiple R-squared:  0.2978,	Adjusted R-squared:  0.2693 
#F-statistic: 10.46 on 3 and 74 DF,  p-value: 8.133e-06

p1 = ggplot(d_params.metadata.sp, 
    aes(
        x = tmax, 
        y = ctmax.ssh, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean maximum temperature of warmest month (",degree~C,")")), 
        y = expression(paste("Maximum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp, 
    aes(
        x = bio9, 
        y = ctmax.ssh, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean temperature of driest quarter (",degree~C,")")), 
        y = expression(paste("Maximum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p2

p3 = ggplot(d_params.metadata.sp, 
    aes(
        x = bio8, 
        y = ctmax.ssh, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean temperature of wettest quarter (",degree~C,")")), 
        y = expression(paste("Maximum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p3

p2 = ggplot(d_params.metadata.sp,
    aes(x = ctmax.ssh)
) +
    geom_histogram(binwidth = 1) +
    labs(
        x = expression(paste("Maximum growth temperature (",degree~C,")"))
    ) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p2

p3 = ggplot(d_params.metadata.sp,
    aes(x = ctmax.boat)
) +
    geom_histogram(binwidth = 0.5) +
    labs(
        x = expression(paste("Boatman maximum growth temperature (",degree~C,")"))
    ) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p3

p4 = ggplot(d_params.metadata.sp, 
    aes(
        x = tmax, 
        y = ctmax.ssh, 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_smooth(method = "lm", alpha = 0.25, color = "black") +
    geom_point(aes(color = state.name)) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_manual(values = c25) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Mean maximum temperature of warmest month (",degree~C,")")), 
        y = expression(paste("Maximum growth temperature (",degree~C,")"))
    ) +
    my_gg_theme
p4

pdf("figures/rTPC_local_adap_comps/ctmax.ssh.pdf", width = 11, height = 5)
p1
p4
p2
p3
dev.off()



###########################################################################
# ## breadth50 ~ Tmax-Tmin (intra-annual), stdev of MAT andor GDD (inter-annual variation)

mod1 = lm(breadth50 ~ I(tmax-tmin)*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)   
#(Intercept)          20.01650    6.49265   3.083  0.00288 **
#I(tmax - tmin)       -0.09336    0.16917  -0.552  0.58269   
#sppNf                -4.93722    6.91480  -0.714  0.47747   
#I(tmax - tmin):sppNf  0.18473    0.18120   1.019  0.31130   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.701 on 74 degrees of freedom
#Multiple R-squared:  0.2142,	Adjusted R-squared:  0.1823 
#F-statistic: 6.722 on 3 and 74 DF,  p-value: 0.0004501

mod2 = lme(fixed = breadth50 ~ I(tmax-tmin)*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same

mod1 = lm(breadth50 ~ HDD4.sd_growing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

#Call:
#lm(formula = breadth50 ~ HDD4.sd_growing * spp, data = d_params.metadata.sp)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-6.5710 -0.6626  0.2771  0.8860  5.0774 
#
#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           16.082107   2.202127   7.303 2.67e-10 ***
#HDD4.sd_growing        0.002797   0.016875   0.166    0.869    
#sppNf                  3.477459   2.451022   1.419    0.160    
#HDD4.sd_growing:sppNf -0.011415   0.018642  -0.612    0.542    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.714 on 74 degrees of freedom
#Multiple R-squared:  0.2029,	Adjusted R-squared:  0.1706 
#F-statistic:  6.28 on 3 and 74 DF,  p-value: 0.0007431

mod2 = lme(fixed = breadth50 ~ HDD4.sd_growing*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same

mod1 = lm(breadth50 ~ HDD4.sd_nongrowing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

colnames(d_params.metadata.sp)

#Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
#(Intercept)              14.50846    2.75195   5.272 1.29e-06 ***
#HDD4.sd_nongrowing        0.06554    0.09228   0.710    0.480    
#sppNf                     4.54921    2.94687   1.544    0.127    
#HDD4.sd_nongrowing:sppNf -0.08557    0.09769  -0.876    0.384    

mod2 = lme(fixed = breadth50 ~ HDD4.sd_nongrowing*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same


mod1 = lm(breadth50 ~ bio2*spp, data = d_params.metadata.sp) # mean diurnal range
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
# ns
mod1 = lm(breadth50 ~ bio3*spp, data = d_params.metadata.sp) #isothermality
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#ns
mod1 = lm(breadth50 ~ bio4*spp, data = d_params.metadata.sp) #T seasonality
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

d_params.metadata.sp$tmax - d_params.metadata.sp$tmin

p1 = ggplot(d_params.metadata.sp, 
    aes(
        y = breadth50, 
        x = tmax-tmin, 
        color = factor(spp, labels = c("N. ditissima", "N. faginata")), 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_point() +
    geom_smooth(method = "lm", alpha = 0.25) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Temperature annual range (",degree~C,")")), 
        y = expression(atop(paste("Niche breadth (",degree~C,")")))
    ) +
    my_gg_theme
p1

p2 = ggplot(d_params.metadata.sp,
    aes(x = breadth50)
) +
    geom_histogram(binwidth = 0.5) +
    labs(
        x = expression(atop(paste("Niche breadth (",degree~C,")")))
    ) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    my_gg_theme
p2

p3 = ggplot(d_params.metadata.sp, 
    aes(
        y = breadth50, 
        x = tmax-tmin, 
        linetype = factor(spp, labels = c("N. ditissima", "N. faginata"))
    )
) +
    geom_smooth(method = "lm", alpha = 0.25, color = "black") +
    geom_point(aes(color = state.name)) +
    scale_linetype_manual(values = c(2,2), guide = "none") +
    scale_color_manual(values = c25) +
    facet_wrap(~factor(spp, labels = c("N. ditissima", "N. faginata")), scale = "free_x") +
    labs(
        x = expression(paste("Temperature annual range (",degree~C,")")), 
        y = expression(atop(paste("Niche breadth (",degree~C,")")))
    ) +
    my_gg_theme
p3


mod1 = lm(breadth50 ~ I(HDD4.mean_growing-HDD4.mean_nongrowing)*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)

pdf("figures/rTPC_local_adap_comps/breadth50.pdf", width = 11, height = 5)
p1
p2
p3
dev.off()


###########################################################################
# ## skewness ~ Tmax

mod1 = lm(skewness ~ tmax*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns

mod2 = lme(fixed = skewness ~ tmax*spp, random = (~1|state.name), data = d_params.metadata.sp)
summary(mod2)
anova(mod2)
#same


###########################################################################
# ## ctmax-topt) ~ Tmax-Tmin (intra-annual), stdev of MAT andor GDD (inter-annual variation)

mod1 = lm(ctmax.ssh-topt ~ I(tmax-tmin)*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#ns

mod1 = lm(ctmax.ssh-topt ~ tmax*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#ns

mod1 = lm(ctmax.ssh-topt ~ HDD4.sd_growing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#ns
#

mod1 = lm(ctmax.ssh-topt ~ HDD4.sd_nongrowing*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#ns
#

mod1 = lm(ctmax.ssh-topt ~ bio8*spp, data = d_params.metadata.sp)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)



#####################
#####################
# we will look at genetic correaltions to topt + t50low (those with sig cors)
# as well as ctmax (both t50low and ctmax have multimodal dists




