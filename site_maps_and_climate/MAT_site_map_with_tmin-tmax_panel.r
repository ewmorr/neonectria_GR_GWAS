library(devtools)
library(reshape2)
library(dplyr)
library(raster)

library(ggplot2)
library(ggmap)

#install_github(repo = "ropensci/prism")
library(prism)

source("library/ggplot_theme.txt")

prism_set_dl_dir("/Users/ericmorrison/PRISM_data")
#options(prism.path = "/Users/ericmorrison/PRISM_data")
#get_prism_normals(type = 'ppt', resolution = '4km', annual = T, keepZip = TRUE)
#get_prism_normals(type = 'tmax', resolution = '4km', annual = T, keepZip = TRUE)
#get_prism_normals(type = 'tmean', resolution = '4km', annual = T, keepZip = TRUE)
#get_prism_normals(type = 'tmin', resolution = '4km', annual = T, keepZip = TRUE)


#Convert raster to point data

new_file<-3#this number corresponds to the row of the file of interest
RS <- pd_stack(prism_archive_ls()[3]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") ##assign projection info

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns


site_coords = read.table("data/sample_metadata/site_coords_for_map.txt", header = T)

minLat=min(site_coords$lat)-2
maxLat=max(site_coords$lat)+2
minLon=min(site_coords$lon)-2
maxLon=max(site_coords$lon)+2


m.df.study_area<-m.df%>%filter(minLat < lat, lat < maxLat, minLon < lon, lon <maxLon)%>%
mutate(MAT = value)%>%
dplyr::select(-value)

min_max = m.df.study_area$MAT %>% range
#min_max[2] = min_max[2]+2
mid_point = (min_max[1]+min_max[2])/2

# need to nudge NH.JG, NH.WF, NH.STR
site_coords

#7    CW1 43.13400 -70.95103      NH.CW  1  0   Nf
#8   DUR1 43.09071 -70.92467      NH.JG  0  1   Nd
#9   DUR2 43.15271 -70.94009      NH.WF  1  1   Nd
#16   STR 43.24958 -71.11645     NH.STR  0  1   Nd
# 
site_coords[7, 3] = -70.825
site_coords[8, 2] = 42.9
site_coords[9, 2] = 43.1
site_coords[9, 3] = -71.05
site_coords[16, 2] = 43.255
site_coords[16, 3] = -71.2

p1 = ggplot()+
    geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=MAT))+
    #geom_point(data=site_coords, aes(x=lon, y = lat, color = state.name), size = 5) +
    #geom_point(data=site_coords, aes(x=lon, y = lat), color = "black", size = 5, shape = 1, show.legend = F) +
    geom_point(data=site_coords, aes(x=lon, y = lat, color = state.name, shape = spp), size = 5) +
    geom_point(
        data=site_coords, 
        aes(x=lon, y = lat), 
        color = "black", size = 5, 
        shape = site_coords$spp.shape.border, 
        show.legend = F
        ) +
    labs(x = "longitude", y = "latitude", shape = "Species sampled", color = "Site") +
    scale_fill_gradient2(expression("Mean annual temperature ("*degree*C*")"),
        low='darkslateblue',
        mid="#d1e5f0",
        high = 'red',
        #mid_point = mid_point
        midpoint=8.5
    ) +
    scale_color_manual(values = c25) +
    scale_shape_manual(values = c(15,17,19), labels = c("both", "N. ditissima", "N. faginata")) +
    guides(
        color = guide_legend(direction = "horizontal", order = 3, title.position = "top"),
        fill = guide_colorbar(direction = "horizontal", order = 1, title.position = "top", barwidth = unit(10, "cm")),
        shape = guide_legend(direction = "horizontal", order = 2, title.position = "top")
    ) +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.77875,0.235),
        legend.key = element_blank(),
        legend.direction = "horizontal"
    )
p1

plot_height = (maxLat-minLat)/2
plot_width = (maxLon-minLon)/2

plot_height.n = plot_height*1.45
plot_width.n = plot_width*1.45

pdf("figures/MAT_site_map.pdf", height = plot_height.n, width = plot_width.n)
print(p1)
dev.off()

###########################
# Climate panel

climate_dat = read.table("data/sample_metadata/sites_climate.txt", header = T)

site.climate = left_join(
    site_coords %>% dplyr::select(Site, state.name), 
    climate_dat %>% dplyr::select(Site, MAT, tmax, tmin, lat)
) %>% tidyr::pivot_longer(
    cols = c(tmin, MAT, tmax),
    names_to = "var",
    values_to = "value"
)

p2 = ggplot(site.climate, 
    aes(x = reorder(state.name, lat), y = value, 
        shape = factor(
            var, 
            levels = c("tmin", "MAT", "tmax"), 
            labels = c("Tmin", "MAT", "Tmax")
        ), 
        fill = state.name,
        color = state.name
    )
) +
    geom_path(aes(group = state.name)) +
    geom_point(size = 4, color = "black") +
    scale_shape_manual(values = c(21,24,22)) +
    scale_fill_manual(values = c25, guide = "none") +
    scale_color_manual(values = c25, guide = "none") +
    labs(y = expression(paste("30-year normals 1981-2010 (", degree~C, ")"))) +
    my_gg_theme +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.045,0.13),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = "#EBEBEB"),
        panel.grid.minor.y = element_line(color = "#F0F0F0")
    )
p2    

library(gridExtra)

pdf("figures/MAT_site_map_w_normals_plot.pdf", height = plot_height.n, width = plot_width.n*1.5)
grid.arrange(p1, p2, ncol = 2, widths = c(0.6667,0.3333))
dev.off()

pdf("figures/MAT_site_map_w_normals_plot.pdf", height = plot_height.n*1.4, width = plot_width.n)
grid.arrange(p1, p2, nrow = 2, heights = c(0.7,0.3))
dev.off()

    
    
    