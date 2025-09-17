library(sn)
library(dplyr)
library(ggplot2)
source("R_scripts/ggplot_theme.txt")
####################################################
# We are currently using a skew-normal distribution
# with no intercept to estimate maximum growth rate,
# optimal growth temp, temperature niche breadth at
# half maximum growth rate, and the lower and upper
# bounds of temperature niche. the distribution 
# was determined after comparison of model fit by
# nls_fit_model_comparisons.R and 
# plots_model_comparisons.R
# GR is originally estimated from growth curves in 
# format_tables_calc_gr.R

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

#look at reps per temp for isolates with higher (bimodal) Topt
odd_ones = c(
    "ADN1 19.1.1",
    "DUR1 1.1.1",
    "GK1j 3.2.1",
    "MEN1 9.1.1",
    "MI1 1.1.1",
    "MI1 2.2.1",
    "MM B2",
    "NJOR2 1.1.1",
    "NJOR2 2.1.2",
    "STR 4.2.2",
    "TSP1 28.2.1"
)
gr_stat.df.metadata.sp %>%
    filter(!is.na(gr.est) ) %>%
    group_by(set, iso_name, temp) %>% 
    summarize(n_reps = n()) %>%
    filter(iso_name %in% odd_ones & temp == 30) %>%
    print(n = Inf)
#all have full replication 

all_temps_names = gr_stat.df.metadata.sp %>%
    filter(!is.na(gr.est) ) %>%
    group_by(iso_name, temp) %>% 
    summarize(n_reps = n()) %>%
    group_by(iso_name) %>%
    summarise(n_temps = n()) %>%
    filter(n_temps > 4) %>%
    pull(iso_name)

gr_stat.df.metadata.sp.full = gr_stat.df.metadata.sp %>%
    filter(iso_name %in% all_temps_names)
nrow(gr_stat.df.metadata.sp)
nrow(gr_stat.df.metadata.sp.full)
################################
#define some starting parameters for a*dsn() model
# based on loess fits and the temp range

# m = Topt (mean x)
# s = (max(x)-min(x))/4
# a = GRmax
# k = 0 #skewness parameter

m = 20
s = (30-10)/4
a = 0.5
k = 0.1

######################
#Functions

# the model fit for skew-normal
f.skew.msak = function(params, x) {
    m <- params[1]; s <- params[2]; a <- params[3]; k = params[4]
    yexp = a*dsn(x = x, xi = m, omega = s, alpha = k) 
    return(yexp)
}

#wrap a formula function and return resid ssq (e.g., for optim)
f.return_ssq = function(params, fctn, x, y){ 
    yexp = fctn(params = params, x = x)
    ssq = sum( (y - yexp)^2 )
    print(paste(params))
    return(ssq)
}

#optim wrapper func
#This is kind of pointless unless we inherit the pars from the environment
# and then run with no parameters, otherwise, same amount of code
f.optim_wrap = function(params, fn, fctn, iters, x, y){
    optim.est = optim(
        par = params,
        fn = fn, #f.return_ssq,
        method = "BFGS",
        control = list(maxit = iters),
        fctn = fctn,
        x = x,
        y = y
    )
    return(optim.est)
}

# niche breadth is the T range at half max GR
# we will calculate yest before the function
# x is the input to the formula function (dense x)
# maxGR estimated by max(yest)
f.breadthHalfMax = function(x, yest) {
    
    #halfMax GR for estimation
    halfMax = max(yest)/2
    
    #break the temp vector and gr vec at the highpoint in GR
    x.low = x[ 1:which.max(yest) ] #lower half
    x.high = x[ which.max(yest):length(x) ] #upper half
    yest.low = yest[ 1:which.max(yest) ] #lower half
    yest.high = yest[ which.max(yest):length(x) ] #upper half
    
    #get the temp at minimum difs between yest and halfmax or max based on absolute values
    Topt = x[ which.max(yest) ]
    Tlow = x.low[ which.min(abs(halfMax - yest.low)) ] 
    Thigh = x.high[ which.min(abs(halfMax - yest.high)) ] 
    
    return(c(Tlow, Topt, Thigh))
}
#End functions
###########################

###########################
# Loop for fits on isolates

#define df
df_len = length(all_temps_names)

fits.df = data.frame(
    iso_name = vector(mode = "character", length = df_len),
    m = vector(mode = "numeric", length = df_len), 
    s = vector(mode = "numeric", length = df_len),
    a = vector(mode = "numeric", length = df_len),# this is mostly useless except for lennon case
    k = vector(mode = "numeric", length = df_len),
    maxGR.obs = vector(mode = "numeric", length = df_len), # avg gr at max temp
    maxGR.exp = vector(mode = "numeric", length = df_len), 
    Tlow = vector(mode = "numeric", length = df_len), 
    Topt = vector(mode = "numeric", length = df_len), 
    Thigh = vector(mode = "numeric", length = df_len), 
    niche_breadth = vector(mode = "numeric", length = df_len),
    
    stringsAsFactors = F
)
#list to store estimates fro each iso
fits_xyest.list = list()


m = 20
s = (30-10)/4
a = 0.5
k = -1

#should maybe try multiple restarts of optim around diff values of k
#and compare ssq of the final pars for different runs to optimize k

#may want to look further into ratkowski

params = c(m,s,a,k)
options(warn=2)
for(iso in 1:length(all_temps_names)){
    
    xy = gr_stat.df.metadata.sp.full %>% filter(iso_name == all_temps_names[iso]) %>% 
        select(temp, gr.est) %>%
        na.omit()
    colnames(xy) = c("x", "y")
    x = xy$x
    y = xy$y
    
    print(iso)
    
    #first try optim
    iters = 100000
    optim.est = optim(
        par = params,
        fn = f.return_ssq,
        method = "BFGS",
        control = list(maxit = iters),#, parscale = c(1,1,1,1)),
        fctn = f.skew.msak,
        x = x,
        y = y
    )
    #check for convergence
    while(optim.est$convergence != 0){
        print("no convergence, retry")
        iters = iters*10
        optim.est = optim(
            par = params,
            fn = f.return_ssq,
            method = "BFGS",
            control = list(maxit = iters),
            fctn = f.skew.msak,
            x = x,
            y = y
        )
    }
    
    #get remaining par estimates
    xest = seq(0,40,0.1)
    yest = f.skew.msak(optim.est$par, x = xest)
    Tniche = f.breadthHalfMax(yest = yest, x = xest)
    #x any ests to list
    fits_xyest.list[[ all_temps_names[iso] ]] = data.frame(x = xest, y = yest)
    
    #assign to df
    fits.df[iso, "iso_name"] = all_temps_names[iso]
    fits.df[iso, c("m", "s", "a", "k")] = optim.est$par
    fits.df[iso, "maxGR.obs"] = data.frame(x,y) %>% 
        group_by(x) %>%
        summarize(mean_y = mean(y) ) %>%
        pull(mean_y) %>%
        max()
    fits.df[iso, "maxGR.exp"] = max(yest)
    fits.df[iso, c("Tlow", "Topt", "Thigh")] = Tniche
    fits.df[iso, "niche_breadth"] = Tniche[3] - Tniche[1]
}

fits.df 
fits_xy.df = bind_rows(fits_xyest.list, .id = "iso_name")

write.table(fits.df, "data/summary_tables/skew_normal_fit_parameters_all_samples.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(fits.df, "data/summary_tables/skew_normal_fits_xy_all_samples.csv", sep = ",", quote = F, col.names = T, row.names = F)

##########################
#PLOT fits versus obs

#join MAT and spp for plotting

#metadata df
spp_MAT = gr_stat.df.metadata.sp.full %>%
    select(iso_name, spp, MAT) %>%
    unique() 
#pull parameters into new df with three points per isolate
#x = tlow, tmean, thigh; y = half max GR, max GR, half max GR

#this rowwise/nest_by operation is a bit more straight forward than 
#pivot_longer select/pull etc given the need to recalc a couple points
pars.list = fits.df %>%
    nest_by(iso_name) %>%
    mutate(
        xy = list(
            data.frame(
                x = c(data$Tlow, data$Topt, data$Thigh) ,
                y = c(data$maxGR.exp/2, data$maxGR.exp, data$maxGR.exp/2) 
            )
        )
    ) 
names(pars.list$xy) =  pars.list$iso_name
pars.df = bind_rows(pars.list$xy, .id = "iso_name")

#join to metadata
pars.df.metadata = left_join(
    pars.df, spp_MAT
)
nrow(pars.df.metadata)
nrow(pars.df)

fits_xy.df.metadata = left_join(
    fits_xy.df, spp_MAT
)
nrow(fits_xy.df.metadata)
nrow(fits_xy.df)

#########################
p = ggplot() +
    geom_point(
        data = pars.df.metadata %>% arrange(rev(spp)) %>% arrange(MAT),
        aes(x = x, y = y),
        color = "blue",
        size = 3,
        alpha = 0.6
    ) +
    geom_line(
        data = fits_xy.df.metadata %>% arrange(rev(spp)) %>% arrange(MAT),
        aes(x = x, y = y),
        color = "black"
    ) +
    geom_point(
        data = gr_stat.df.metadata.sp.full %>% arrange(rev(spp)) %>% arrange(MAT),
        aes(y = gr.est, x = temp, shape = as.factor(rep), color = spp ),
        position = position_jitter(width = 0.2, height = 0)
    ) +
    facet_wrap(spp~iso_name) +
    labs(
        y = expression(paste("Growth rate (mm day"^-1, ")")), 
        x = expression(paste("Growth temperature (",degree,C, ")")),
        shape = "Replicate",
        color = "Species"
    ) +
    scale_color_manual(values = cbPalette) +
    theme_bw() 

pdf("figures/gr_fits_skew_normal.pdf", width = 12, height = 12)
p
dev.off()

##############################################################################
# parameter m is not estimating Topt correctly (i.e., temperature at maxGR)
# for highly skewed data. Should pull directly from estimated GR~T curves

# Also need to increase T range for estimation (from original 3-35C)
# there is one sample (VAS1 7.1) where Thigh does not fall within 35C max
# and is therefore incorrect (bc our function pull closest difference)

# Above issues are now addressed within this script

