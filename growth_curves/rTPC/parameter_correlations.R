library(dplyr)
library(GGally)
library(ggplot2)
source("library/ggplot_theme.txt")
library(vegan)
library(ggrepel)


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


p <- ggpairs(
  d_params.metadata.sp %>% filter(spp == "Nf") %>% select(rmax, topt, ctmax.ssh, ctmax.boat, skewness, t50low, t50high, breadth50),
  lower = list(continuous = wrap("points", alpha = 0.6, size = 1.5)),
  upper = list(continuous = "cor", combo = "box_no_facet", discrete = "count", na = "na"),
  diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag"),
)
p

colnames(d_params.metadata.sp)
pairs(d_params.metadata.sp %>% select(rmax, topt, ctmax.ssh, ctmax.boat, skewness, t50low, t50high, breadth50))

pdf("figures/rTPC/parameter_corrs.pdf", width = 12, height = 10)
p
dev.off()


# pca

# Run PCA (center + scale recommended)
pca <- prcomp(d_params.metadata.sp %>% filter(spp == "Nf") %>% select(rmax, topt, ctmax.ssh, ctmax.boat, skewness, t50low, t50high, breadth50), center = TRUE, scale. = TRUE)

# Scores (samples in PC space)
scores <- as.data.frame(pca$x)
scores.meta = cbind(d_params.metadata.sp %>% filter(spp == "Nf") %>% select(iso_name), scores)
write.csv(scores.meta, "data/rTPC/growth_params_PCA.csv", quote = F, row.names = F)

# Loadings (variable contributions to PCs)
loadings <- as.data.frame(pca$rotation)
loadings$var <- rownames(loadings)

# Scale loadings for plotting as arrows
loadings_plot <- loadings
loadings_plot$PC1 <- loadings$PC1 * max(abs(scores$PC1))
loadings_plot$PC2 <- loadings$PC2 * max(abs(scores$PC2))

# Percent variance explained
var_expl <- round(100 * summary(pca)$importance[2, 1:2], 1)

pdf("figures/rTPC/parameter_PCA.scree.pdf")
plot(round(100 * summary(pca)$importance[2,],1))
dev.off()

# PCA biplot
p1 = ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, size = 2) +  # sample scores
  geom_segment(data = loadings_plot,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "dark grey") +
  geom_text_repel(data = loadings_plot,
                  aes(x = PC1*1.1, y = PC2*1.1, label = var),
                  color = "#3c3c3c",
                  size = 5,
                  segment.color = NA) +   # no leader lines
  my_gg_theme +
  labs(
       x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)")
    )
p1

pdf("figures/rTPC/parameter_PCA.pdf", width = 10, height = 8)
p1
dev.off()


######################################
#####climate cors with PCs
######################################

d_params.metadata.sp.Nf = left_join(d_params.metadata.sp %>% filter(spp == "Nf"), scores.meta)


library(nlme)

#PC1
#
mod1 = lme(fixed = PC1 ~ MAT, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns

mod1 = lme(fixed = PC1 ~ HDD4.mean_growing, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns

mod1 = lme(fixed = PC1 ~ HDD4.mean_nongrowing, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns

mod1 = lme(fixed = PC1 ~ length.mean_growing, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns

mod1 = lme(fixed = PC1 ~ tmin, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns


#PC2

mod1 = lme(fixed = PC2 ~ tmax, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns

mod1 = lme(fixed = PC2 ~ breadth50, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns

mod1 = lme(fixed = PC2 ~ HDD4.sd_growing, random = (~1|state.name), data = d_params.metadata.sp.Nf)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
anova(mod1)
#ns


