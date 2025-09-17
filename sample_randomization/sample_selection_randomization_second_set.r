#Randomize samples within site to a max of four samples per site
require(dplyr)

#setwd("~/GARNAS_neonectria_GR_GWAS/")

sample_df = read.table("data/sequenced_samples_for_growth_rate_randomization_NF_06092022.passed_qual.txt", header = T, sep = "\t")

#sample(0:x, 4, replace = F) #Random selection on array

sample_df.site_levels = levels(as.factor(sample_df$state.region))

sample.list = list()

for(i in 1:length(sample_df.site_levels)){
    
    temp.df = sample_df %>% filter(state.region == sample_df.site_levels[i])
    n_samps = nrow(temp.df)
    if(n_samps > 4){
        rand_indices = sample(1:n_samps, 4, replace = F)
        sample.list[[sample_df.site_levels[i]]] = temp.df[rand_indices,]
    }else{
        sample.list[[sample_df.site_levels[i]]] = temp.df
    }
}

#do.call("rbind", sample.list) #apparently this is much slower than bind_rows(), time to update

sample.random_selection = bind_rows(sample.list)

write.table(sample.random_selection, file = "data/sequenced_samples_random_select_NF_06092022.txt", sep = "\t", quote = F, row.names = F)
