#  Repo for analysis of growth rate data using GWAS approach

# Data dir is "~/GARNAS_neonectria_GR_GWAS"


```
format_tables_calc_gr.r # processing raw growth rate data and calculation of growth rate
rTPC/rTPC.multiple_curves.r # calculation and comparison of model fits
rTPC/spp_comparisons.r # comparing modeled growth parameters between spp
rTPC/local_adaptation_lms.r # comparison of growth metrics to climate data
```

compare standard deviation within pops to annual sd
compare topt to tmean, winter GDD and summer GDD
compare ctmax to tmax and summer GDD; and ctmin to tmin and winter GDD
ANOVA of pops on GRmax, Topt (and other metrics)

## spp comparisons
- sig difs between Nf and Nd in:
    - t50high (high temp at half max GR), Nf higher
    - GRmax, Nd higher
    - CTmax, Nf higher
    - Niche breadth (at half max GR), Nf higher

## local adaptation
- test for an effect of species; in all cases where there was a significant relationship between growth metrics and climate there was a significant interaction with species such that Nf was significant and Nd was not
- sig correlation between:
    - Topt and MAT
    - Topt and growing season GDD
    - t50low and nongrowing season GDD
    - ctmax appears to have multimodal dist in Nf

### Test for genetic signal in Topt, t50low, and ctmax with LFMM
