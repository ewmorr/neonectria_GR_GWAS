library(lfmm)
library(dplyr)
library(ggplot2)
library(scales)
source("library/ggplot_theme.txt")

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

# join to sample level metadata. Just the species and sequence label
# sample metadata
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[3:5] = c("iso_name", "sequence_id", "spp")
d_params.metadata.sp = left_join(
    d_params.metadata,
    sample_ID_map %>% select(iso_name, spp, sequence_id) %>% unique,
    by = "iso_name"
)
head(d_params.metadata.sp)
sum(complete.cases(d_params.metadata.sp)) == nrow(d_params.metadata.sp)
d_params.metadata.sp.Nf = d_params.metadata.sp %>% filter(spp == "Nf")
rownames(d_params.metadata.sp.Nf) = d_params.metadata.sp.Nf$sequence_id
d_params.metadata.sp.Nf
d_params.metadata.sp.Nf %>% filter(iso_name == "ANF1 10.2")

d_params.metadata.sp.Nf %>% pull(Site) %>% unique()


########################################################################
########################Completed data set up###########################
########################################################################


#PED sample IDs
sample_ids = read.table("data/genotype_dat/FINAL_snp.mac_ge2.biallele.retainNG152.gwas_analyses.sampleIDs", header = F)
sum(rownames(d_params.metadata.sp.Nf) %in% sample_ids$V1) == nrow(d_params.metadata.sp.Nf)
sample_ids.Nf = sample_ids$V1[sample_ids$V1 %in% d_params.metadata.sp.Nf$sequence_id]
length(sample_ids.Nf) == nrow(d_params.metadata.sp.Nf)
# sort
d_params.metadata.sp.Nf.sorted = d_params.metadata.sp.Nf[sample_ids.Nf,]
sum(complete.cases(d_params.metadata.sp.Nf.sorted)) == nrow(d_params.metadata.sp.Nf.sorted)

row_ids = which(sample_ids$V1 %in% d_params.metadata.sp.Nf.sorted$sequence_id)

#The genotype data can simply be read in as a matrix (according the docs)
#OR can try loading LEA and using readLfmm()
Y = as.matrix(read.table("data/genotype_dat/FINAL_snp.mac_ge2.biallele.retainNG152.gwas_analyses.lfmm", header = F))
Y.filtered = Y[row_ids,]
nrow(Y.filtered)
ncol(Y.filtered)

SNP_pos = read.table("data/genotype_dat/FINAL_snp.mac_ge2.biallele.retainNG152.gwas_analyses.recode.map")
nrow(SNP_pos)
SNP_pos = SNP_pos[c(1,4)]
colnames(SNP_pos) = c("scaffold", "position")


###################
###################
###################
#MAC filter
ref_sum = colSums(Y.filtered == 0)
alt_sum = colSums(Y.filtered == 1)

#filter for MAC ge 3
nrow(Y.filtered)
minMAC = 3

sum(ref_sum < minMAC)
# 101504
sum(alt_sum < minMAC)
ncol(Y.filtered)
# 424676
ncol(Y.filtered) - sum(ref_sum < minMAC)
# 269836

which(ref_sum < minMAC) 
which(alt_sum < minMAC)# there are none
rm_cols = which(ref_sum < minMAC)

length(rm_cols)#154840
ncol(Y.filtered)
length(rm_cols)/ncol(Y.filtered)
# 0.3646074
ncol(Y.filtered)-length(rm_cols)
# 269836
Y.filteredMAC = Y.filtered[,-rm_cols]
SNP_pos.filteredMAC = SNP_pos[-rm_cols,]
ncol(Y.filteredMAC)
nrow(SNP_pos.filteredMAC)

#filter for greater than 25% missing data
#sum(colSums(Y.filteredMAC == 9)/nrow(Y.filteredMAC) > 0.25)
#rm_cols = which(colSums(Y.filteredMAC == 9)/nrow(Y.filteredMAC) > 0.25)
#Y.filteredNA = Y.filteredMAC[,-rm_cols]

plot(colSums(Y == 9)/nrow(Y))
plot(colSums(Y == 0))
plot(colSums(Y == 1))
###################
###################
###################

Y.filtered = Y.filteredMAC
SNP_pos = SNP_pos.filteredMAC

ncol(Y.filtered)
nrow(SNP_pos)

#
#principal components analysis for K
pc = prcomp(Y.filtered)
plot((pc$sdev^2)/sum(pc$sdev^2), xlab = 'PC', ylab = "% variance explained")
(pc$sdev^2)/sum(pc$sdev^2)
points(4,pc$sdev[7]^2, type = "h", lwd = 3, col = "blue")
str(pc)
plot(pc$x[,2] ~ pc$x[,1])
plot(pc$x[,3] ~ pc$x[,1])
plot(pc$x[,3] ~ pc$x[,2])
plot(pc$x[,4] ~ pc$x[,1])

#Y = Y.filteredNA


#######################
#######################
#topt
#######################
#######################

#variable for test
X = (d_params.metadata.sp.Nf.sorted$topt)
#X = (sample_metadata.site_info[,39:43])


class(X)
class(Y.filtered)
sd(X)
mean(X)
#LFMM ridge

mod.lfmm = lfmm_ridge(Y = Y.filtered, X = X, K = 3, lambda = 1) #using K = 3 based on PCA 
str(mod.lfmm)

pv <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm,
calibrate = "gif")
str(pv)

#Example plots
plot(-log10(pv$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "SNP", ylab = "-Log P",
col = "grey")

plot((pv$B),
pch = 19,
cex = .3,
xlab = "SNP", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv$score^2, na.rm = T)/0.456
lambda #0.9951334
adj.p.values = pchisq(pv$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv$calibrated.pvalue) #this looks conservative
hist(pv$pvalue)
#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/0.9, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 0.9
adj.p.values = pchisq(pv$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)
length(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/genotype_dat/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

#Join with actual positiion and chromosome
nrow(SNP_pos)
ncol(Y.filtered)
pv.with_pos = data.frame(calibrated.p = pv$calibrated.pvalue, effect_size = pv$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.with_pos )$calibrated.p, 0.025, na.rm = T) #removed the filter by 100000, can do this later for plotting but it does not seem good to do before outlier ID; %>% filter(length > 100000)
# make an outlier column in the data.frame
pv.with_pos <- pv.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
#note that once we filtered for MAC we have no more NA p values
pv.with_pos %>% group_by(outlier) %>% tally()
pv.with_pos %>% filter(is.na(outlier)) %>% head
#6746
#FDR correction
#This is based on the auto calibartion
pv.with_pos$FDR.p = p.adjust(pv.with_pos$calibrated.p, method = "fdr", n = length(pv.with_pos$calibrated.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig) %>% tally()

#1  SNPs identified as significant after FDR correction


#FDR correction
#This is based on the manual GIF adjustment
pv.with_pos$FDR.p.man = p.adjust(pv.with_pos$man.adj.p, method = "fdr", n = length(pv.with_pos$man.adj.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig.man) %>% tally()
#4 (5 at 0.85, 4 at 0.9)

pv.topt.with_pos = pv.with_pos



pv.topt.with_pos %>% head()
pv.topt.with_pos %>% filter(outlier == "outlier") %>% nrow()
#6746 
# but this is not the same thing as sig these are 2.5% outliers

####################
#ggplots

#Basic plot
ggplot(
    pv.topt.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8)
    #axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#Colored by outliers
ggplot(
    pv.topt.with_pos %>% filter(length > 100000 ), 
    aes(x = position/10^6, y = calibrated.p, color = outlier)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black")) +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom"
)

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.topt.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.topt.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/LFMM/topt.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/topt.man_cal.png", width = 1080, height = 240)
p2
dev.off()

#######################
#######################
#t50low
#######################
#######################

#variable for test
X = d_params.metadata.sp.Nf.sorted$t50low

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 3)
mod.lfmm.t50low = lfmm_ridge(Y = Y.filtered, X = X, K = 3)

pv.t50low <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.t50low,
calibrate = "gif")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.t50low$score^2, na.rm = T)/0.456
lambda #1.048667
adj.p.values = pchisq(pv.t50low$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.t50low$calibrated.pvalue)
hist(pv.t50low$pvalue)
#IN THIS CASE THE CALCULATED VALUES look conservative

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.t50low$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.t50low$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.t50low$score^2/1.0, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS conservative
#However, the lower values have a correct distribution under null model
#Try at GIF = 1 (both 0.95 and 1 have the same after FDR)
adj.p.values = pchisq(pv.t50low$score^2/1, df = 1, lower = FALSE)
hist(adj.p.values)

nrow(SNP_pos)
nrow(scf_lens)
#Join with actual positiion and chromosome
pv.t50low.with_pos = data.frame(calibrated.p = pv.t50low$calibrated.pvalue, effect_size = pv.t50low$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.t50low.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.t50low.with_pos <- pv.t50low.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.t50low.with_pos %>% group_by(outlier) %>% tally()

#6746

#FDR correction
pv.t50low.with_pos = pv.t50low.with_pos %>% filter(!is.na(calibrated.p))

pv.t50low.with_pos$FDR.p = p.adjust(pv.t50low.with_pos$calibrated.p, method = "fdr", n = length(pv.t50low.with_pos$calibrated.p))
range(pv.t50low.with_pos$FDR.p, na.rm = T)
pv.t50low.with_pos <- pv.t50low.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.t50low.with_pos %>% group_by(FDR.sig) %>% tally()

#0 at 0.05

#FDR correction
#This is based on the manual GIF adjustment
pv.t50low.with_pos$FDR.p.man = p.adjust(pv.t50low.with_pos$man.adj.p, method = "fdr", n = length(pv.t50low.with_pos$man.adj.p))
pv.t50low.with_pos <- pv.t50low.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.t50low.with_pos %>% group_by(FDR.sig.man) %>% tally()
#0 at 1, 4 at 0.95

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.t50low.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.t50low.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/LFMM/t50low.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/t50low.man_cal.png", width = 1080, height = 240)
p2
dev.off()

#######################
#######################
#ctmax
#######################
#######################

#variable for test
X = d_params.metadata.sp.Nf.sorted$ctmax.ssh

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.ctmax = lfmm_ridge(Y = Y.filtered, X = X, K = 3)

pv.ctmax <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.ctmax,
calibrate = "gif")


#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.ctmax$score^2, na.rm = T)/0.456
lambda #1.265398
adj.p.values = pchisq(pv.ctmax$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.ctmax$calibrated.pvalue)
hist(pv.ctmax$pvalue)
#IN THIS CASE THE CALCULATED VALUES lookmaybe slightly conservative but pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ctmax$score^2/1.35, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ctmax$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ctmax$score^2/1.2, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good, but more conservative

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.2
adj.p.values = pchisq(pv.ctmax$score^2/1.2, df = 1, lower = FALSE)
hist(adj.p.values)


#Join with actual positiion and chromosome
nrow(SNP_pos)
nrow(pv.ctmax$calibrated.pvalue)
pv.ctmax.with_pos = data.frame(calibrated.p = pv.ctmax$calibrated.pvalue, effect_size = pv.ctmax$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.ctmax.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.ctmax.with_pos <- pv.ctmax.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.ctmax.with_pos %>% group_by(outlier) %>% tally()


#FDR correction
pv.ctmax.with_pos$FDR.p = p.adjust(pv.ctmax.with_pos$calibrated.p, method = "fdr", n = length(pv.ctmax.with_pos$calibrated.p))
pv.ctmax.with_pos <- pv.ctmax.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.ctmax.with_pos %>% group_by(FDR.sig) %>% tally()

#5

#FDR correction manual adjustment (lambda = 1)
pv.ctmax.with_pos$FDR.p.man = p.adjust(pv.ctmax.with_pos$man.adj.p, method = "fdr", n = length(pv.ctmax.with_pos$man.adj.p))
pv.ctmax.with_pos <- pv.ctmax.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.ctmax.with_pos %>% group_by(FDR.sig.man) %>% tally()

#6 at 1.2, 7 at 1.15

####################
#ggplots


#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.ctmax.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8)
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()

)
p1

p2 = ggplot(
    pv.ctmax.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8)
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()
)
p2

png("figures/LFMM/ctmax.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/ctmax.man_cal.png", width = 1080, height = 240)
p2
dev.off()


#######################
#######################
#skewness
#######################
#######################

#variable for test
X = d_params.metadata.sp.Nf.sorted$skewness

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 3)
mod.lfmm.skewness = lfmm_ridge(Y = Y.filtered, X = X, K = 3)

pv.skewness <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.skewness,
calibrate = "gif")

#Example plots
plot(-log10(pv.skewness$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "-Log P",
col = "grey")

plot((pv.skewness$B),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.skewness$score^2, na.rm = T)/0.456
lambda #1.03922
adj.p.values = pchisq(pv.skewness$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.skewness$calibrated.pvalue)
hist(pv.skewness$pvalue)
#IN THIS CASE THE CALCULATED VALUES look pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.skewness$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.skewness$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.skewness$score^2/1.0, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS conservative
#However, the lower values have a correct distribution under null model
#Try at GIF = 1 (both 0.95 and 1 have the same after FDR)
adj.p.values = pchisq(pv.skewness$score^2/1, df = 1, lower = FALSE)
hist(adj.p.values)

nrow(SNP_pos)
nrow(scf_lens)
#Join with actual positiion and chromosome
pv.skewness.with_pos = data.frame(calibrated.p = pv.skewness$calibrated.pvalue, effect_size = pv.skewness$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.skewness.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.skewness.with_pos <- pv.skewness.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.skewness.with_pos %>% group_by(outlier) %>% tally()

#6746

#FDR correction
pv.skewness.with_pos = pv.skewness.with_pos %>% filter(!is.na(calibrated.p))

pv.skewness.with_pos$FDR.p = p.adjust(pv.skewness.with_pos$calibrated.p, method = "fdr", n = length(pv.skewness.with_pos$calibrated.p))
range(pv.skewness.with_pos$FDR.p, na.rm = T)
pv.skewness.with_pos <- pv.skewness.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.skewness.with_pos %>% group_by(FDR.sig) %>% tally()

#1606 at 0.05

#FDR correction
#This is based on the manual GIF adjustment
pv.skewness.with_pos$FDR.p.man = p.adjust(pv.skewness.with_pos$man.adj.p, method = "fdr", n = length(pv.skewness.with_pos$man.adj.p))
pv.skewness.with_pos <- pv.skewness.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.skewness.with_pos %>% group_by(FDR.sig.man) %>% tally()
#1771 at 1, 2097 at 0.95

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.skewness.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.skewness.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/LFMM/skewness.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/skewness.man_cal.png", width = 1080, height = 240)
p2
dev.off()



#######################
#######################
# PC1 (growth parametr PCA)
#######################
#######################

pc_dat = read.csv("data/rTPC/growth_params_PCA.csv")
pc_dat

d_params.metadata.sp.Nf.sorted.pcs = left_join(d_params.metadata.sp.Nf.sorted, pc_dat)
#variableread.csv2(#variable for test
X = d_params.metadata.sp.Nf.sorted.pcs$PC1

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 3)
mod.lfmm.grpc1 = lfmm_ridge(Y = Y.filtered, X = X, K = 3)

pv.grpc1 <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.grpc1,
calibrate = "gif")

#Example plots
plot(-log10(pv.grpc1$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "-Log P",
col = "grey")

plot((pv.grpc1$B),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.grpc1$score^2, na.rm = T)/0.456
lambda #1.263475
adj.p.values = pchisq(pv.grpc1$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.grpc1$calibrated.pvalue)
hist(pv.grpc1$pvalue)
#IN THIS CASE THE CALCULATED VALUES look pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.grpc1$score^2/1.35, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.grpc1$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.grpc1$score^2/1.2, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS conservative
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.2 
adj.p.values = pchisq(pv.grpc1$score^2/1.2, df = 1, lower = FALSE)
hist(adj.p.values)

nrow(SNP_pos)
nrow(scf_lens)
#Join with actual positiion and chromosome
pv.grpc1.with_pos = data.frame(calibrated.p = pv.grpc1$calibrated.pvalue, effect_size = pv.grpc1$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.grpc1.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.grpc1.with_pos <- pv.grpc1.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.grpc1.with_pos %>% group_by(outlier) %>% tally()

#6746

#FDR correction
pv.grpc1.with_pos = pv.grpc1.with_pos %>% filter(!is.na(calibrated.p))

pv.grpc1.with_pos$FDR.p = p.adjust(pv.grpc1.with_pos$calibrated.p, method = "fdr", n = length(pv.grpc1.with_pos$calibrated.p))
range(pv.grpc1.with_pos$FDR.p, na.rm = T)
pv.grpc1.with_pos <- pv.grpc1.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.grpc1.with_pos %>% group_by(FDR.sig) %>% tally()

#6 at 0.05

#FDR correction
#This is based on the manual GIF adjustment
pv.grpc1.with_pos$FDR.p.man = p.adjust(pv.grpc1.with_pos$man.adj.p, method = "fdr", n = length(pv.grpc1.with_pos$man.adj.p))
pv.grpc1.with_pos <- pv.grpc1.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.grpc1.with_pos %>% group_by(FDR.sig.man) %>% tally()
#19 at 1.2, 25 at 1.15

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.grpc1.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.grpc1.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/LFMM/grpc1.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/grpc1.man_cal.png", width = 1080, height = 240)
p2
dev.off()


#######################
#######################
# PC2 (growth parametr PCA)
#######################
#######################

pc_dat = read.csv("data/rTPC/growth_params_PCA.csv")
pc_dat

d_params.metadata.sp.Nf.sorted.pcs = left_join(d_params.metadata.sp.Nf.sorted, pc_dat)
#variableread.csv2(#variable for test
X = d_params.metadata.sp.Nf.sorted.pcs$PC2

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 3)
mod.lfmm.grpc2 = lfmm_ridge(Y = Y.filtered, X = X, K = 3)

pv.grpc2 <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.grpc2,
calibrate = "gif")


#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.grpc2$score^2, na.rm = T)/0.456
lambda #0.9671996
adj.p.values = pchisq(pv.grpc2$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.grpc2$calibrated.pvalue)
hist(pv.grpc2$pvalue)
#IN THIS CASE THE CALCULATED VALUES look pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.grpc2$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.grpc2$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.grpc2$score^2/0.9, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS conservative
#However, the lower values have a correct distribution under null model
#Try at GIF = 0.9
adj.p.values = pchisq(pv.grpc2$score^2/0.9, df = 1, lower = FALSE)
hist(adj.p.values)

nrow(SNP_pos)
nrow(scf_lens)
#Join with actual positiion and chromosome
pv.grpc2.with_pos = data.frame(calibrated.p = pv.grpc2$calibrated.pvalue, effect_size = pv.grpc2$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.grpc2.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.grpc2.with_pos <- pv.grpc2.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.grpc2.with_pos %>% group_by(outlier) %>% tally()

#6615

#FDR correction
pv.grpc2.with_pos = pv.grpc2.with_pos %>% filter(!is.na(calibrated.p))

pv.grpc2.with_pos$FDR.p = p.adjust(pv.grpc2.with_pos$calibrated.p, method = "fdr", n = length(pv.grpc2.with_pos$calibrated.p))
range(pv.grpc2.with_pos$FDR.p, na.rm = T)
pv.grpc2.with_pos <- pv.grpc2.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.grpc2.with_pos %>% group_by(FDR.sig) %>% tally()

#16 at 0.05

#FDR correction
#This is based on the manual GIF adjustment
pv.grpc2.with_pos$FDR.p.man = p.adjust(pv.grpc2.with_pos$man.adj.p, method = "fdr", n = length(pv.grpc2.with_pos$man.adj.p))
pv.grpc2.with_pos <- pv.grpc2.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.grpc2.with_pos %>% group_by(FDR.sig.man) %>% tally()
#34 at 0.9, 62 at 0.85

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.grpc2.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.grpc2.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/LFMM/grpc2.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/grpc2.man_cal.png", width = 1080, height = 240)
p2
dev.off()


#######################
#######################
# MAT (climate)
#######################
#######################


#variableread.csv2(#variable for test
X = d_params.metadata.sp.Nf.sorted$MAT

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 3)
mod.lfmm.mat = lfmm_ridge(Y = Y.filtered, X = X, K = 3)

pv.mat <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.mat,
calibrate = "gif")


#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.mat$score^2, na.rm = T)/0.456
lambda #1.183144
adj.p.values = pchisq(pv.mat$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.mat$calibrated.pvalue)
hist(pv.mat$pvalue)
#IN THIS CASE THE CALCULATED VALUES look pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.mat$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.mat$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.mat$score^2/1.1, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS conservative
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.15 
adj.p.values = pchisq(pv.mat$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

nrow(SNP_pos)
nrow(scf_lens)
#Join with actual positiion and chromosome
pv.mat.with_pos = data.frame(calibrated.p = pv.mat$calibrated.pvalue, effect_size = pv.mat$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.mat.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.mat.with_pos <- pv.mat.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.mat.with_pos %>% group_by(outlier) %>% tally()

#6746

#FDR correction
pv.mat.with_pos = pv.mat.with_pos %>% filter(!is.na(calibrated.p))

pv.mat.with_pos$FDR.p = p.adjust(pv.mat.with_pos$calibrated.p, method = "fdr", n = length(pv.mat.with_pos$calibrated.p))
range(pv.mat.with_pos$FDR.p, na.rm = T)
pv.mat.with_pos <- pv.mat.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.mat.with_pos %>% group_by(FDR.sig) %>% tally()

#256 at 0.05

#FDR correction
#This is based on the manual GIF adjustment
pv.mat.with_pos$FDR.p.man = p.adjust(pv.mat.with_pos$man.adj.p, method = "fdr", n = length(pv.mat.with_pos$man.adj.p))
pv.mat.with_pos <- pv.mat.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.mat.with_pos %>% group_by(FDR.sig.man) %>% tally()
#261 at 1.15

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.mat.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.mat.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/LFMM/mat.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/mat.man_cal.png", width = 1080, height = 240)
p2
dev.off()


#######################
#######################
# HDD4.growing (climate)
#######################
#######################


#variableread.csv2(#variable for test
X = d_params.metadata.sp.Nf.sorted$HDD4.mean_nongrowing

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 3)
mod.lfmm.hdd4g = lfmm_ridge(Y = Y.filtered, X = X, K = 3)

pv.hdd4g <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.hdd4g,
calibrate = "gif")


#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.hdd4g$score^2, na.rm = T)/0.456
lambda #1.297482
adj.p.values = pchisq(pv.hdd4g$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.hdd4g$calibrated.pvalue)
hist(pv.hdd4g$pvalue)
#IN THIS CASE THE CALCULATED VALUES look pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.hdd4g$score^2/1.35, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.hdd4g$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.hdd4g$score^2/1.2, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS conservative
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.15 
adj.p.values = pchisq(pv.hdd4g$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values)

nrow(SNP_pos)
nrow(scf_lens)
#Join with actual positiion and chromosome
pv.hdd4g.with_pos = data.frame(calibrated.p = pv.hdd4g$calibrated.pvalue, effect_size = pv.hdd4g$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.hdd4g.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.hdd4g.with_pos <- pv.hdd4g.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.hdd4g.with_pos %>% group_by(outlier) %>% tally()

#6746

#FDR correction
pv.hdd4g.with_pos = pv.hdd4g.with_pos %>% filter(!is.na(calibrated.p))

pv.hdd4g.with_pos$FDR.p = p.adjust(pv.hdd4g.with_pos$calibrated.p, method = "fdr", n = length(pv.hdd4g.with_pos$calibrated.p))
range(pv.hdd4g.with_pos$FDR.p, na.rm = T)
pv.hdd4g.with_pos <- pv.hdd4g.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.hdd4g.with_pos %>% group_by(FDR.sig) %>% tally()

#2327 at 0.05

#FDR correction
#This is based on the manual GIF adjustment
pv.hdd4g.with_pos$FDR.p.man = p.adjust(pv.hdd4g.with_pos$man.adj.p, method = "fdr", n = length(pv.hdd4g.with_pos$man.adj.p))
pv.hdd4g.with_pos <- pv.hdd4g.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.hdd4g.with_pos %>% group_by(FDR.sig.man) %>% tally()
#2563 at 1.15

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.hdd4g.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig, alpha = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.hdd4g.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/LFMM/hdd4_growing.auto_cal.png", width = 1080, height = 240)
p1
dev.off()

png("figures/LFMM/hdd4_growing.man_cal.png", width = 1080, height = 240)
p2
dev.off()



#######################
#######################
#######################
#######################
#######################
#######################
#write tables

pv.topt.with_pos %>% head
pv.t50low.with_pos %>% head
pv.ctmax.with_pos %>% head

write.table(pv.topt.with_pos, "data/LFMM/topt.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.t50low.with_pos, "data/LFMM/t50low.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.ctmax.with_pos, "data/LFMM/ctmax.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.grpc1.with_pos, "data/LFMM/grpc1.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.grpc2.with_pos, "data/LFMM/grpc2.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.mat.with_pos, "data/LFMM/mat.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.hdd4g.with_pos, "data/LFMM/hdd4_growing.lfmm.txt", quote = F, row.names = F, sep = "\t")

#Nice aligned plot of all three
library(gtable)
library(gridExtra)
library(grid)

############################################
##READ THE ABOVE TABLES BACK IN TO RUN BELOW

p1 = ggplot(
    pv.topt.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    #axis.title.y = element_blank(),
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)

p2 = ggplot(
    pv.t50low.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)

p3 = ggplot(
    pv.ctmax.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    #axis.title.y = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8)
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()
)


########################
#manual CALIPBRATED P-VALUES
plots = list(p1, p2, p3)
plots = list(p1, p3)

grobs <- lapply(plots, ggplotGrob)
# for gridExtra < v2.3, use do.call(gridExtra::rbind.gtable, grobs)
# for gridExtra >= v2.3 use:
g <- do.call(gridExtra::gtable_rbind, grobs)

#set relative panel heights
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(rep(1,2), "null")


grid.newpage()
grid.draw(g)


png("figures/LFMM/topt-t50low-ctmax_manual_GIF.png", width = 1080, height = 420)
grid.newpage()
grid.draw(g)
dev.off()

png("figures/LFMM/topt-ctmax_manual_GIF.png", width = 1080, height = 320)
grid.newpage()
grid.draw(g)
dev.off()


###############################################################################
###############################################################################
###############################################################################

# PC plots

p1 = ggplot(
    pv.grpc1.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    #axis.title.y = element_blank(),
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)


p2 = ggplot(
    pv.grpc2.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man, alpha = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
scale_alpha_manual(values = c(0.5, 1), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("-log"[10],"(",italic(P),")"))) +
theme(
    #axis.title.y = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8)
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()
)


########################
#manual CALIPBRATED P-VALUES
plots = list(p1, p2)

grobs <- lapply(plots, ggplotGrob)
# for gridExtra < v2.3, use do.call(gridExtra::rbind.gtable, grobs)
# for gridExtra >= v2.3 use:
g <- do.call(gridExtra::gtable_rbind, grobs)

#set relative panel heights
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(rep(1,2), "null")


grid.newpage()
grid.draw(g)


png("figures/LFMM/grpc1-grpc2_manual_GIF.png", width = 1080, height = 320)
grid.newpage()
grid.draw(g)
dev.off()

