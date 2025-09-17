library(sn)
library(dplyr)

###################################################################
# This script implements comparisons of four models for GR ~ T
# Gaussian according to Lennon et al. 2012, Gaussian with intercept
# skew normal, skew normal plus intercept. See nls_fit_tests.R for
# model formulation. 20 isolates with all growth curves having 
# r^2 > 0.98 selected randomly for model comparisons. Parameter
# estimates of GRmax, Topt, and skew compared between models as well
# as GRmax estimate - observed max GR for each indvidual isolate.
# Ideal model would have *consistent* deviation from observed max GR
# and also have consistent relationship to parameter estimates of 
# other models.

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

####################################
#Select a set of isolates for comps 

#Filter out growth curves with >0.97-0.99 r^2
#count number remaining reps per isolate/temp
#filter out levels with less than 3 reps except 5 and 35 C where 2 reps



gr_stat.df.metadata.sp %>% 
    filter(gr.rsq > 0.98) %>%
    group_by(set, iso_name, temp) %>%
    summarize(reps = sum(!is.na(gr.est))) %>%
    filter(
        (temp %in% c(10,15,20,25,30) & reps >= 3) |
        (temp %in% c(5, 35) & reps >= 1) 
    ) %>% 
    group_by(set, iso_name) %>%
    summarize(temps = n()) %>%
    filter(temps > 6) %>%
    print(n = Inf)

#0.97 R2 gives 47 isos, 0.98 34, 0.99 8
#Use 0.98, excluding those in third filter. (why were these filtered before? (!iso_name %in% c("DUR2 1.2", "ITH1 1.2", "MEN2 1.2", "MEN2 6.4") ))
# All curves of GR over T look reasonably tight
# i.e., small standard error, and there are no obvious outliers
#

names_good = gr_stat.df.metadata.sp %>% 
    filter(gr.rsq > 0.98) %>%
    group_by(set, iso_name, temp) %>%
    summarize(reps = sum(!is.na(gr.est))) %>%
    filter(
        (temp %in% c(10,15,20,25,30) & reps >= 3) |
        (temp %in% c(5, 35) & reps >= 1) 
    ) %>% 
    group_by(set, iso_name) %>%
    summarize(temps = n()) %>%
    filter(temps > 4) %>%
    pull(iso_name)

set.seed(12345)
isos_sample = sample(names_good, 25, replace = F)

################################
#define some starting parameters 
# based on loess fits and the temp range

# m = Topt (mean x)
# s = (max(x)-min(x))/4
# a = GRmax
# k = 0 #skewness parameter
# b = min(y)

#note this order is retained for all models
# but not all params go to every model
m = 20
s = (30-10)/4
a = 0.5
k = 0
b = 0.1

######################
# Functions for fits

#model names
fs = c(
    "f.gaus.msa",
    "f.gausInt.msab",
    "f.skew.msak",
    "f.skewInt.msakb"
)

#number of params for calculating AIC
#this is number of estimated parameters plus 1 (error)
fs.n_params = list(4,5,5,6)
names(fs.n_params) = fs    

#list for model functions
fs.list = list()

#lennon gaussian
fs.list[[fs[1]]] = f.gaus.msa = function(params, x)  { 
    m <- params[1]; s <- params[2]; a <- params[3]; k = params[4]; b <- params[5];
    yexp = a*exp(-abs((x-m)/s)^2)
    return(yexp)
}

#Gaussian with intercept
fs.list[[fs[2]]] = f.gausInt.msab = function(params, x) {
    m <- params[1]; s <- params[2]; a <- params[3]; k = params[4]; b <- params[5];
    yexp = a*exp(-0.5*((x-m)/s)^2) + b
    return(yexp)
}

#skewed normal 
fs.list[[fs[3]]] = f.skew.msak = function(params, x) {
    m <- params[1]; s <- params[2]; a <- params[3]; k = params[4]; b <- params[5];
    yexp = a*dsn(x = x, xi = m, omega = s, alpha = k) 
    return(yexp)
}

#skewed normal with intercept
fs.list[[fs[4]]] = f.skewInt.msakb = function(params, x) {
    m <- params[1]; s <- params[2]; a <- params[3]; k = params[4]; b <- params[5];
    yexp = a*dsn(x = x, xi = m, omega = s, alpha = k) + b
    return(yexp)
}

fs.list 

#wrap a formula function and return ssq (e.g., for optim)
#need function input arg different than optim fn arg
f.return_ssq = function(params, fctn, x, y){ 
    yexp = fctn(params = params, x = x)
    ssq = sum( (y - yexp)^2 )
    return(ssq)
}

#AIC calculation
f.aic = function(n_obs, n_params, rss){
    temp.mle = rss/(n_obs - n_params)
    aic = 2*n_params + n_obs*log(temp.mle)
    return(aic)
}

#End fit functions
###################

#################
#optim test

#params = c(m,s,a,k,b)
#y = gr_stat.df.metadata.sp %>% filter(iso_name == isos_sample[1]) %>% pull(gr.est) %>% na.omit
#x = gr_stat.df.metadata.sp %>% filter(iso_name == isos_sample[1]) %>% pull(temp) %>% na.omit

#f.return_ssq(params = params, fctn = f.gaus.msa, x = x, y = y)

#optim(
#    par = params,
#    fn = f.return_ssq,
#    method = "BFGS",
#    control = list(maxit = 10000),
#    fctn = fs.list[[fs[1]]],
#    x = x,
#    y = y
#)

# will want to condition on $convergence
# set maxit in object and increase if $convergence != 0

# for loop can just input all pars for the wrapped formulas. 
# The ones not used are ignored (i.e., don't affect ssq)
# need to declare all pars in formula fn though

# End test
###################

###################
#optim wrapper func
#This is kind of pointless unless we inherit the pars from the environment
# and then run with no parameters, otherwise, same amount of code
f.optim_wrap = function(params, fn, fctn, iters, x, y){
    optim.est = optim(
        par = params,
        fn = fn, #f.return_ssq,
        method = "BFGS",
        control = list(maxit = iters),
        fctn = fctn, #fs.list[[f.name]],
        x = x,
        y = y
    )
    return(optim.est)
}
###################

###################
#functions for parameter estimates

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
    
    #get the temp at minimum difs between yest and halfmax based on absolute values
    Tlow = x.low[ which.min(abs(halfMax - yest.low)) ] 
    Thigh = x.high[ which.min(abs(halfMax - yest.high)) ] 
    
    return(c(Tlow, Thigh))
}

###########################
# Loop for fits on isolates

#define df
df_len = length(isos_sample)*length(fs)

fits.df = data.frame(
    iso_name = vector(mode = "character", length = df_len),
    model_name = vector(mode = "character", length = df_len),
    a = vector(mode = "numeric", length = df_len), # this is also mostly useless except for lennon case
    m = vector(mode = "numeric", length = df_len),
    s = vector(mode = "numeric", length = df_len),
    k = vector(mode = "numeric", length = df_len),
    b = vector(mode = "numeric", length = df_len), # this probably isn't important but let's retain
    maxGR.obs = vector(mode = "numeric", length = df_len), # avg gr at max temp
    maxGR.exp = vector(mode = "numeric", length = df_len), 
    T.low = vector(mode = "numeric", length = df_len), 
    T.high = vector(mode = "numeric", length = df_len), 
    niche_breadth = vector(mode = "numeric", length = df_len),
    rss = vector(mode = "numeric", length = df_len),
    aic = vector(mode = "numeric", length = df_len),
    
    stringsAsFactors = F
)

#loop
counter = 0
params = c(m,s,a,k,b)

for(iso in 1:length(isos_sample)){
    
    xy = gr_stat.df.metadata.sp %>% filter(iso_name == isos_sample[iso]) %>% 
        select(temp, gr.est) %>%
        na.omit()
    colnames(xy) = c("x", "y")
    x = xy$x
    y = xy$y

    print(iso)
    
    for(mod in 1:length(fs)){
        f.name = fs[mod]
        iters = 10000
        counter = counter + 1
        
        #first try
        optim.est = optim(
            par = params,
            fn = f.return_ssq,
            method = "BFGS",
            control = list(maxit = iters),
            fctn = fs.list[[f.name]],
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
                fctn = fs.list[[f.name]],
                x = x,
                y = y
            )
        }
        
        #get remaining par estimates
        xest = seq(5,35,0.1)
        yest = fs.list[[f.name]](optim.est$par, x = xest)
        T.niche = f.breadthHalfMax(yest = yest, x = xest)
        rss = f.return_ssq(params = optim.est$par, fctn = fs.list[[f.name]], x = x, y = y)
        aic = f.aic(n_obs = length(y), n_params = fs.n_params[[f.name]], rss = rss)
        
        #assign to df
        fits.df[counter, "iso_name"] = isos_sample[iso]
        fits.df[counter, "model_name"] = f.name
        fits.df[counter, c("m", "s", "a", "k", "b")] = optim.est$par
        fits.df[counter, "maxGR.obs"] = data.frame(x,y) %>% 
            group_by(x) %>%
            summarize(mean_y = mean(y) ) %>%
            pull(mean_y) %>%
            max()
        fits.df[counter, "maxGR.exp"] = max(yest)
        fits.df[counter, c("T.low", "T.high")] = T.niche
        fits.df[counter, "niche_breadth"] = T.niche[2] - T.niche[1]
        fits.df[counter, "rss"] = rss
        fits.df[counter, "aic"] = aic
    }
}

fits.df 

fits.df$model_name = gsub("^f\\.", "", fits.df$model_name) %>% 
    gsub("\\.[msakb]+", "", x = ., perl = T)

write.table(fits.df, "data/summary_tables/model_fits_comps.csv", sep = ",", quote = F, col.names = T, row.names = F)
