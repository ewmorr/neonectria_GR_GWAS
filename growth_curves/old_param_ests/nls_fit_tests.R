library(sn)
library(dplyr)

source("growth_curves/make_site_metadata.r")
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
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[c(3, 5)] = c("iso_name", "spp")
gr_stat.df.metadata.sp = left_join(
    gr_stat.df.metadata,
    sample_ID_map %>% select(iso_name, spp) %>% unique,
    by = "iso_name"
)
nrow(gr_stat.df.metadata.sp)

######################
#Maynard is using nls and dsn to estimate a "standard 3-parameter skew normal distribution

#dsn has parameters
# x a vetor of quantiles
# p vector of probabilities
#xi vector of location parameters
#omega vector of scale parameters
# etc

# Lennon et al. 2012 provides a formla for the reponse where tau is the shape of the curve, so dsn in maynards application

# gr ~ GRmax * ( exp( -abs( (temp - Topt)/sigma  ) )^tau   )

# parameters to be estimated
#GRmax (max growth rate)
#Topt (optimal temperature for growth)
#sigma (the rate that respiration declines as a strain moves away from Topt)

# niche breadth (b) is then defined as
# b = sigma( -log10(x)^(1/tau) )
# in Lennon and Maynard x = 0.5, or 50% of maximum gr

#Need to then give initial points for GRmax, Topt, sigma

#estimate parameters from loess smooths of gr ~ temp
GRmax = 0.5
Topt = 22.5
sd = 5 #range in temp/4 

################
# nls formula can be either input directly or as a function. function is probably cleaner
# normal distribution formula from https://stackoverflow.com/questions/61150560/nonlinear-least-squares-regression-of-skewed-normal-distribution-in-r-or-any-la


############################
#Starting with normal distribution

#model of normal distribution 
f <- function(x, theta)  { 
    m <- theta[1]; s <- theta[2]; a <- theta[3]; b <- theta[4];
    a*exp(-0.5*((x-m)/s)^2) + b
}

#version for optim
f.opt <- function(x, y, theta)  { 
    m <- theta[1]; s <- theta[2]; a <- theta[3]; b <- theta[4];
    yexp = a*exp(-0.5*((x-m)/s)^2) + b
    ssq = sum((y - yexp)^2)
    return(ssq)
}

# m = mean of x (Topt)
# s = standard deviation
# a = peak y (GRmax) *AFTER subtracting b (so need to add b for real value) and should estimate as max(y)-min(y)
# b (y intercept)

#get a representative iso_name for testing
test_name = "ASH2 3.1.1"
iso_quants = gr_stat.df.metadata.sp %>% filter(iso_name == test_name) %>% select(gr.est) %>% quantile(na.rm = T) %>% unname

y = gr_stat.df.metadata.sp %>% filter(iso_name == test_name) %>% pull(gr.est)
y
x = gr_stat.df.metadata.sp %>% filter(iso_name == test_name) %>% pull(temp)
x
plot(y ~ x)

#estimate starting values
Topt = x[which.max(y)]
GRmax = max(y)-min(y)
sd = (max(x)-min(x))/4
b = min(y)

fit <- nls(
    y ~ f(x,c(m,s,a,b)), 
    data.frame(x,y), 
    start = list(
        m = Topt, 
        s = sd, 
        a = GRmax, 
        b = b
    ),
    na.omit
)

summary(fit)$parameters["m", 1:2] #21.7
summary(fit)$parameters["a", 1:2] #0.27
summary(fit)$parameters["b", 1:2] #0.21
summary(fit)$parameters["s", 1:2] #4

#plot the fit
par(mfrow=c(1,1))
plot(c(x,0),c(y,f(coef(fit)["m"],coef(fit))), main="Data", type="n",
     xlab="GR", ylab="Temp")
curve(f(x, coef(fit)), add=TRUE, col="Red", lwd=2)
points(x,y, pch=19)

#Now with optim
gaus_fit.optim = optim(
    par = c(Topt, sd, GRmax, b),
    fn = f.opt,
    x = x,
    y = y
)

gaus_fit.optim$par
summary(fit)$parameters
#They are the same. Should use optim bc nls fails in the dsn case

######################
#Define formula based on Lennon
f <- function(x, theta)  { 
    m <- theta[1]; s <- theta[2]; a <- theta[3];
    a*exp(-abs((x-m)/s)^2)
}

# m = mean of x (Topt)
# s = standard deviation
# a = peak y (GRmax); because we do not include a y-intercept term the model estimates the real max
#note that in Lennon the ^2 is denoted as tau and modifies the shape of the function

fit <- nls(
    y ~ f(x,c(m,s,a)), 
    data.frame(x,y), 
    start = list(
        m = Topt, 
        s = sd, 
        a = GRmax
    ),
    na.omit
)

summary(fit)$parameters["m", 1:2] #20.85
summary(fit)$parameters["a", 1:2] # 4386
#summary(fit)$parameters["b", 1:2]
summary(fit)$parameters["s", 1:2] # 11.7

#plot the fit
par(mfrow=c(1,1))
plot(c(x,0),c(y,f(coef(fit)["m"],coef(fit))), main="Data", type="n",
     xlab="Temp", ylab="GR")
curve(f(x, coef(fit)), add=TRUE, col="Red", lwd=2)
points(x,y, pch=19)

###################
#From here Lennon estimates niche breadth (b) based on 
# b = sd * ( -log10(x) )^(1/tau) (where tau equals 2 in the case of normal dist
# x is a modifier that gives the values of x at that proportion of GRmax (or proportion of a)

b = summary(fit)$parameters["s", 1] * (-log10(0.5) )^(1/2)
b
b = summary(fit)$parameters["s", 1] * (-log10(0.00000000000000000001) )^(1/2)
b

#
m=Topt
s=5
a=GRmax
f(x,c(m,s,a))

#######################
# Skewed normal based on https://stackoverflow.com/questions/61150560/nonlinear-least-squares-regression-of-skewed-normal-distribution-in-r-or-any-la

f <- function(x, theta)  { 
    m <- theta[1]; s <- theta[2]; a <- theta[3]; b <- theta[4]; k <- theta[5]
    a*exp(k*((x - m))/s - sqrt(((x - m))/s*((x - m))/s+1) ) + b
}

a = GRmax
m = Topt
s = (max(x)-min(x))/4
b = min(y)
k = 0 #skewness parameter

plot(f(x, c(m, s, a, k, b)) ~ x)

fit <- nls(y ~ f(x,c(m,s,a,b, k)), data.frame(x,y), start=list(m=m, s=s, a=a, b=b, k=k))

summary(fit)$parameters["m", 1:2] #24
summary(fit)$parameters["a", 1:2] #0.8
summary(fit)$parameters["b", 1:2] #0.2
summary(fit)$parameters["s", 1:2] #2.15
summary(fit)$parameters["k", 1:2] #-0.5

#plot the fit
par(mfrow=c(1,1))
plot(c(x,0),c(y,f(coef(fit)["m"],coef(fit))), main="Data", type="n",
     xlab="Temp", ylab="GR", ylim = c(0,0.6))
curve(f(x, coef(fit)), add=TRUE, col="Red", lwd=2)
points(x,y, pch=19)

f(x, coef(fit))
#Can find x vals for given y by increasing the density of x
test_temps = seq(5,35,0.1)

f(seq(5,35,0.1) , coef(fit)) %>% max


  halfGRmax = (summary(fit)$parameters["a", 1] - summary(fit)$parameters["b", 1])/2

f(seq(5,35,0.011) , coef(fit) ) %>% max

####################
#Had something set incorrectly yesterday. The above is now running *but* GRmax seems high ... note that GRmax is now a - b not a + b as in previous formula
summary(fit)$parameters["a", 1] - summary(fit)$parameters["b", 1]

#actually a is not giving the max based on calculating the real values
f(seq(5,35,0.011) , coef(fit) ) %>% max
max.ind = f(seq(5,35,0.011) , coef(fit) ) %>% which.max()
( seq(5,35,0.011) )[max.ind]


###################
#Trying with dsn instead of self defined formula
#using optim as indicated in Maynard

a = GRmax
m = Topt
s = (max(x)-min(x))/4
b = min(y)
k = 0 #skewness parameter

#need to wrap dsn in a new function in order to pass paramters 
#return SSQ from dsn to obs

plot(dsn(x = x, xi = m, omega = s, alpha = k) ~ x)

f.opt = function(param,x,y){
    skewPts = param[5]*dsn(x = x, xi = param[1], omega = param[2], alpha = param[3]) 
    return(sum((y-skewPts)^2))
}

f.opt(c(m, s, k,b,a),x,y)


output.optim = optim(
    par = c(m, s, k, b,a),
    fn = f.opt,
    method = "BFGS",
    control = list(maxit = 10000),
    x = x,
    y = y
)
output.optim
#m = 19.9 no b or a #with a and no b this becomes 20.59 #with a and b becomes 21.5
#s = 2.7
#k = 60.6
#b = 0.21
#a = 2.7

pts.dsn = output.optim$par[5]*dsn(x = seq(0,40,0.1), xi = output.optim$par[1], omega = output.optim$par[2], alpha = output.optim$par[3]) 

ggplot() +
    geom_line(data = data.frame(x = seq(0,40,0.1), y = pts.dsn), aes(x = x, y = y)) +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    theme_bw()


#############
#Trying the above fit with nls
a = GRmax
m = Topt
s = (max(x)-min(x))/4
b = min(y)
k = 0.1 #skewness parameter

fit = nls(
    y ~ a*dsn(x = x, xi = m, omega = s, alpha = k)+b,
    start = list(a = a, m = m, s = s, b = b, k = k),
    data = data.frame(x = x, y = y),
    na.action = na.omit
)

#The below gives a normal distribution with mean (xi) 0 and sd (omega) 1 by setting alpha (the skew paramter) to 0
plot(dsn(x = seq(-3,3, 0.1), xi=0, omega=1, alpha=0, log=FALSE) ~ seq(-3,3, 0.1))
#negative alpha gives left skew
plot(dsn(x = seq(-3,3, 0.1), xi=0, omega=1, alpha= -1, log=FALSE) ~ seq(-3,3, 0.1))
#positive alpha gives right skew
plot(dsn(x = seq(-3,3, 0.1), xi=0, omega=1, alpha= 1, log=FALSE) ~ seq(-3,3, 0.1))


plot(2*(dsn(x = seq(-3,3, 0.1), xi=0, omega=1, alpha=0, log=FALSE) ) ~ seq(-3,3, 0.1) )
     
#The GR data is right-skewed

#selm from the sn package
mod = selm(
    formula = y ~ x,
    family = "SN",
    data = data.frame(x = x, y = y)
)
plot(mod)
summary(mod)
fitted(mod)


ggplot() +
    geom_line(data = data.frame(x = x, y = fitted(mod)), aes(x = x, y = y)) +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    theme_bw()
