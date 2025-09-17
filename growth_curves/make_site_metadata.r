library(dplyr)

#fam_info = read.table("data/sample_metadata/out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
#sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
#colnames(sample_ID_map)[3:4] = c("iso_name", "sample")
site_info = read.table("data/sample_metadata/site_coords.txt", header = T) # 
site_climate = read.table("data/sample_metadata/sites_climate.txt", header = T) # need to join by state.name because some sites are grouped (e.g., MEN1 and MEN2 are both the same location/coords)
site_GDD = read.table("data/sample_metadata/site_climate.GDD.txt", header = T)
worldclim = read.csv("data/sample_metadata/worldclim.csv")

#
site_GDD = left_join(site_GDD, site_info %>% select(Site, state.name)) #for this table not all of the subsites are included. need to map all of the existing sites to state.name and then join to the larger table by state.name while excluding site. This gives the subsites the same climate info as the main sites which is appropriate

site_GDD.state = left_join(
    site_GDD,
    site_info %>% select(Site, state.name)
) %>% left_join(.,
    worldclim %>% select(-lat, -lon)
) %>% filter(!Site %in% c("LWF1", "ADS2", "ADN2", "MEN2", "MM"))

# The filtered sites are duplicate values 

#best option may be to make site to state.name table and then join by site
site_metadata = left_join(
    site_info %>% select(Site, state.name) %>% filter(!Site %in% c("LWF1", "ADS2", "ADN2", "MEN2", "MM")),
    site_climate %>% select(Site, ppt, tmax, MAT, tmin) %>% filter(!Site %in% c("LWF1", "ADS2", "ADN2", "MEN2", "MM"))
) %>% left_join(
    .,
    site_GDD.state %>% select(-Site)
)
site_metadata

#Now join to growth rate df by site to avoid multiple matches to state.name
