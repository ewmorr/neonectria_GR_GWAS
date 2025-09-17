library(dplyr)
library(ggplot2)
source("R_scripts/ggplot_theme.txt")

fits.df = read.csv("data/summary_tables/model_fits_comps.csv")
colnames(fits.df)

#######################
#Compare predicted v obs max GR across models
p1 = ggplot(fits.df, aes(x = maxGR.exp, y = maxGR.obs)) +
    geom_point() +
    facet_wrap(~model_name) +
    geom_smooth(method = "lm", se = F) +
    geom_abline(slope = 1) +
    labs(x = "Predicted max growth rate", y = "Observed max growth rate") +
    my_gg_theme
p1
#intercept models consistently over estimate at high end
#should prob compare slopes and r2

#fun newish dplyr way to run function on subset with nest_by and mutate
lm_models = fits.df %>%
    nest_by(model_name) %>%
    mutate(mod = list(lm(maxGR.obs ~ maxGR.exp, data = data))) #need list() to transform lm obj to a vector
lm_models.coefs = lm_models %>% 
    summarise(
        rsq = round(summary(mod)$r.squared, 2), 
        intercept = round(coef(mod)[1], 2),
        slope = round(coef(mod)[2], 2)
    )
lm_models.labels = lm_models.coefs %>%
    mutate(
        rsq.lab = sprintf("italic(R^2) == %.2f", rsq),
        formula = sprintf("italic(y) == %.2f * italic(x) %+.2f ", slope, intercept)
    )


#Compare predicted v obs max GR across models
p1 = ggplot(fits.df, aes(x = maxGR.exp, y = maxGR.obs)) +
    geom_point() +
    geom_text(
        data = lm_models.labels, 
        x = 0.45, 
        y = 0.4, 
        parse = TRUE, 
        hjust = 0,
        aes(label = formula )
    ) +
    geom_text(
        data = lm_models.labels, 
        x = 0.45, 
        y = 0.39, 
        parse = TRUE, 
        hjust = 0,
        aes(label = rsq.lab )
    ) +
    facet_wrap(~model_name) +
    geom_smooth(method = "lm", se = F) +
    geom_abline(slope = 1) +
    labs(x = "Predicted max growth rate", y = "Observed max growth rate") +
    my_gg_theme
p1

p2 = ggplot(fits.df, aes(x = model_name, y = rss)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.1, height = 0)) +
    labs(x = "Model", y = expression(italic(RSS)) )  +
    my_gg_theme
p2

#comparing aic
p3 = ggplot(fits.df, aes(x = model_name, y = aic)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.1, height = 0)) +
    labs(x = "Model", y = "AIC" )  +
    my_gg_theme
p3

mod.aic = aov(aic ~ model_name, data = fits.df)
summary(mod.aic)
#all models are basically equivalent in terms of AIC
#note there is a group of isolates with higher aic for all models (not anymore)
# this is the set with fewer obs (most likely)

pdf("figures/GR_estimate_model_comparisons.pdf")
p1
p2
p3
dev.off()
