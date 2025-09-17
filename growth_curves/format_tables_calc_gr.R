library(dplyr)
library(lubridate)
library(ggplot2)

data_path = "data/processed_data"
files_vec = list.files(path=data_path, pattern = "*.csv", full.names = T)
length(files_vec)

set_rep_temp.list = basename(files_vec) %>% 
    sub("set", "", x = .) %>%
    sub("rep", "", x = .) %>%
    sub("C", "", x = .) %>%
    sub(".csv", "", x = .) %>%
    strsplit("_")

files_list = vector(mode = "list", length = length(files_vec))

for(i in 1: length(files_vec)){
    files_list[[i]] = read.csv(files_vec[i], header = F)
}

files_list
# times are imported as military
#note that there are some NA cols at the end of some tables due to conversion from the excel data sheets (these will need to be removed obv)

#also importing data from repeated growth curves
redos = read.csv("data/processed_data/redos/repeats.csv", header = F)
redos.rep_temp = read.csv("data/processed_data/redos/repeats.rep_temp.csv", header = F)

redosTwo = read.csv("data/processed_data/redos/repeats.09032024.csv", header = F)
redosTwo.rep_temp = read.csv("data/processed_data/redos/repeats.09032024.rep_temp.csv", header = F)

##################################
#testing for date time processing
#and indexing growth rate data

#dt = files_list[[1]][1:3,-1] %>% t() %>% data.frame %>% na.omit

#the na.omit removes na rows so subsequent operations on the full df should be done on the nrow of this df

#nrow(dt)
#dt
#colnames(dt) = files_list[[1]][1:3,1]

#start and end times

#mdy_hm(paste(dt$Date, dt$`Time Start`, " ") ) 
#mdy_hm(paste(dt$Date, dt$`Time End`, " ") ) 
#dt$start = mdy_hm(paste(dt$Date, dt$`Time Start`, " ") ) 
#dt$end = mdy_hm(paste(dt$Date, dt$`Time End`, " ") ) 
#length(dt$start)
#diff(dt$start) #difftime object
#int_diff(dt$start) %>% int_length # length of interval in seconds; numeric
#int_diff(dt$start) %>% int_length / (3600*24) #days

#append zero for starting measurement and cumsum hours

#dt$start.hr_sum = c(0, int_diff(dt$start) %>% int_length / (3600*24)) %>% cumsum
#dt$end.hr_sum = c(0, int_diff(dt$end) %>% int_length / (3600*24)) %>% cumsum
#dt

#mean of measurement time

#( dt$start.hr_sum + dt$end.hr_sum ) / 2
#dt$mean.hr_sum = ( dt$start.hr_sum + dt$end.hr_sum ) / 2

####
# testing gr data

#i = 1
#number_measures = nrow(dt) #number of measurement times across the dataset 
#iso_name = files_list[[i]][ 4:nrow(files_list[[i]]) , 1 ] %>% t() %>% c()

#gr = files_list[[i]][ 4:nrow(files_list[[i]]) , 2:(number_measures+1) ] %>% #add one to account for excluding names column 
#    t() %>%
#    apply(., 2, as.numeric) %>% #need to convert from character string
#    data.frame() #convert to df bc there are some tables with funny values that muck up matrix-based
#gr
#is.numeric(gr)
#colnames(gr) = iso_name
#gr
#object.size(gr)

#END testing
###############

#####################
#data processing loop

#lists for storing processed data
dt.list = vector(mode = "list", length = length(files_vec))
gr.list = vector(mode = "list", length = length(files_vec))

#when added set 1 rep 3 file 14 started throwing a NAs introduced by coercion error
#set_rep_temp.list[[14]]
#set 1, rep 3, T25
#one cell had comma instead of period. If this error crops up again check out the csv for commas

options(warn = 2) #set warnings to interrupt (like an error)
for(i in 1: length(files_list)){
    
    ###############
    #datetime data
    
    #the na.omit removes na rows. Subsequent operations on the full df should be done on the length of this df
    dt = files_list[[i]][1:3,-1] %>% t() %>% data.frame %>% na.omit
    number_measures = nrow(dt)
    colnames(dt) = files_list[[i]][1:3,1]
    
    #set start and end of measurement window datetime
    dt$start = mdy_hm(paste(dt$Date, dt$`Time Start`, " ") ) 
    dt$end = mdy_hm(paste(dt$Date, dt$`Time End`, " ") ) 
    
    #calculate intervals in numeric hours
    #append zero for starting measurement and cumsum days intervals
    dt$start.day_sum = c(0, int_diff(dt$start) %>% int_length / (3600*24)) %>% cumsum
    dt$end.day_sum = c(0, int_diff(dt$end) %>% int_length / (3600*24)) %>% cumsum
    
    #mean of measurement window
    dt$mean.day_sum = ( dt$start.day_sum + dt$end.day_sum ) / 2
    
    #save to list 
    dt.list[[i]] = dt
    
    ##################
    #gr data
    
    #number_measures = nrow(dt) #number of measurement times across the dataset
    iso_name = files_list[[i]][ 4:nrow(files_list[[i]]) , 1 ] %>% t() %>% c()
    gr = files_list[[i]][ 4:nrow(files_list[[i]]) , 2:(number_measures+1) ] %>% #add one to account for excluding names column 
        t() %>%
        apply(., 2, as.numeric) %>%
        data.frame()
    #print(paste(i, is.numeric(gr)) )
    colnames(gr) = iso_name
    gr = gr[colnames(gr) != ""] #again, some blank cols import from excel
    gr.list[[i]] = gr
}
i #if throw an error check file
files_vec[i] #if throw an error check file
options(warn=0)
#End main data processing loop
#########################

dt.list
gr.list
set_rep_temp.list

########################
#repeats growth curves processing loop
#for now we dump these into three dfs (dt, gr, rep_temp)

{
    ###########
    #process dt
    redos.dt = redos[1:3,-1] %>% t() %>% data.frame %>% na.omit
    colnames(redos.dt) = redos[1:3,1]
    #set start and end of measurement window datetime
    redos.dt$start = mdy_hm(paste(redos.dt$Date, redos.dt$`Time Start`, " ") ) 
    redos.dt$end = mdy_hm(paste(redos.dt$Date, redos.dt$`Time End`, " ") ) 
    
    #calculate intervals in numeric hours
    #append zero for starting measurement and cumsum days intervals
    redos.dt$start.day_sum = c(0, int_diff(redos.dt$start) %>% int_length / (3600*24)) %>% cumsum
    redos.dt$end.day_sum = c(0, int_diff(redos.dt$end) %>% int_length / (3600*24)) %>% cumsum
    
    #mean of measurement window
    redos.dt$mean.day_sum = ( redos.dt$start.day_sum + redos.dt$end.day_sum ) / 2
    
    ###########
    #process gr
    redos.gr = redos[ 4:nrow(redos) , 2:ncol(redos) ] %>%
        t() %>%
        apply(., 2, as.numeric) %>%
        data.frame()
    redos.iso_name = redos[4:nrow(redos), 1]
    
    nrow(redos.rep_temp)
    ncol(redos.gr)
    
    ###########
    #Need to create a temp list subset by rep and temp to collect all the isos
    # that have the same rep and temp name
    
    redos.rows.list = list()
    unique.sets = unique(redos.rep_temp[c("V1", "V2", "V3")])
    for(i in 1:nrow(unique.sets)){
        redos.rows.list[[ paste(unique.sets[i,], collapse = "_") ]] = 
            which(
                apply(
                    redos.rep_temp, 1, 
                    function(x) all(x == unique.sets[i,])
                    )
                )
    }
    
    # run through the list of unique rep_temps and add each set of growth  
    # curves to the list. These will be sets of growth curves added to the 
    # main lists by set, rep, temp 
    
    main.len = length(files_vec) #length of primary data lists
    for(i in 1:length(redos.rows.list)){
        
        set_rep_temp.list[[main.len + i]] = 
            unlist(strsplit(names(redos.rows.list)[[i]], "_") ) 
        
        
        gr.list[[main.len + i]] = setNames(
            data.frame(redos.gr[,redos.rows.list[[i]]] ), 
            redos.iso_name[redos.rows.list[[i]]]
        )
        
        dt.list[[main.len + i]] = redos.dt #all get the same dt
    }
}
set_rep_temp.list

########################
#repeats growth curves processing loop TAKE TWO
#for now we dump these into three dfs (dt, gr, rep_temp)

{
    ###########
    #process dt
    redosTwo.dt = redosTwo[1:3,-1] %>% t() %>% data.frame %>% na.omit
    colnames(redosTwo.dt) = redosTwo[1:3,1]
    #set start and end of measurement window datetime
    redosTwo.dt$start = mdy_hm(paste(redosTwo.dt$Date, redosTwo.dt$`Time Start`, " ") ) 
    redosTwo.dt$end = mdy_hm(paste(redosTwo.dt$Date, redosTwo.dt$`Time End`, " ") ) 
    
    #calculate intervals in numeric hours
    #append zero for starting measurement and cumsum days intervals
    redosTwo.dt$start.day_sum = c(0, int_diff(redosTwo.dt$start) %>% int_length / (3600*24)) %>% cumsum
    redosTwo.dt$end.day_sum = c(0, int_diff(redosTwo.dt$end) %>% int_length / (3600*24)) %>% cumsum
    
    #mean of measurement window
    redosTwo.dt$mean.day_sum = ( redosTwo.dt$start.day_sum + redosTwo.dt$end.day_sum ) / 2
    
    ###########
    #process gr
    redosTwo.gr = redosTwo[ 4:nrow(redosTwo) , 2:ncol(redosTwo) ] %>%
        t() %>%
        apply(., 2, as.numeric) %>%
        data.frame()
    redosTwo.iso_name = redosTwo[4:nrow(redosTwo), 1]
    
    nrow(redosTwo.rep_temp)
    ncol(redosTwo.gr)
    
    ###########
    #Need to create a temp list subset by rep and temp to collect all the isos
    # that have the same rep and temp name
    
    redosTwo.rows.list = list()
    unique.sets = unique(redosTwo.rep_temp[c("V1", "V2", "V3")])
    for(i in 1:nrow(unique.sets)){
        redosTwo.rows.list[[ paste(unique.sets[i,], collapse = "_") ]] = 
            which(
                apply(
                    redosTwo.rep_temp, 1, 
                    function(x) all(x == unique.sets[i,])
                    )
                )
    }
    
    # run through the list of unique rep_temps and add each set of growth  
    # curves to the list. These will be sets of growth curves added to the 
    # main lists by set, rep, temp 
    
    main.len = length(files_vec) #length of primary data lists
    for(i in 1:length(redosTwo.rows.list)){
        
        set_rep_temp.list[[main.len + i]] = 
            unlist(strsplit(names(redosTwo.rows.list)[[i]], "_") ) 
        
        
        gr.list[[main.len + i]] = setNames(
            data.frame(redosTwo.gr[,redosTwo.rows.list[[i]]] ), 
            redosTwo.iso_name[redosTwo.rows.list[[i]]]
        )
        
        dt.list[[main.len + i]] = redosTwo.dt #all get the same dt
    }
}
set_rep_temp.list
##########################
#Stats and plots loop

#Can index plot list by iso name, then by rep.temp to print in order
#Just start with blank list

#table will need cols for iso, set, rep, temp, x_est, x_p, x_rsq, and n_obs. Table length is ncol gr.list[[i]] (should be consistent across list)
lapply(gr.list, ncol)
#set2 has dif number of cols bc dif number of tables
#it's OK bc setting length by ncol in the loop

#make table in top level of loop and populate set rep temp info at once then proceed to iso level of loop

#set up empty lists
gr_stat.list = vector(mode = "list", length = length(files_vec))
plots.list = list()

for(i in 1: length(gr.list) ){
    iso_set = set_rep_temp.list[[i]][1]
    iso_rep = set_rep_temp.list[[i]][2]
    iso_temp = set_rep_temp.list[[i]][3]
    
    #data.frame for storing gr stats
    gr.df = data.frame(
        iso_name = vector(mode = "character", ncol(gr.list[[i]])),
        set = rep(iso_set, ncol(gr.list[[i]])),
        rep = rep(iso_rep, ncol(gr.list[[i]])),
        temp = rep(iso_temp, ncol(gr.list[[i]])),
        gr.est = vector(mode = "numeric", ncol(gr.list[[i]])),
        gr.p = vector(mode = "numeric", ncol(gr.list[[i]])),
        gr.rsq = vector(mode = "numeric", ncol(gr.list[[i]])),
        n_obs = vector(mode = "numeric", ncol(gr.list[[i]])),
        stringsAsFactors = F
    )
    
    #days as x
    x = dt.list[[i]]$mean.day_sum
    
    for(j in 1: ncol(gr.list[[i]])){
        iso_name = colnames(gr.list[[i]])[j]
        y = gr.list[[i]][,j] * 7.5 #mult by 1.5 for conversion to mm (from tube grads)
        #THE ABOVE NEEDS TO CHANGE. EACH 0.2 ML MARKER IS 1.5 MM NOT A FULL 1 ML GRAD
        #INSTEAD DO THE GRAD MEASUREMENTS TIMES 5*1.5(7.5) I.E. THERE ARE 7.5 ML PER 1 GRAD
        #DOUBLE CHECKED BOTH THE 0.2 ML AND 1 ML CONVERSION ON 2/6
        
        #A few isos have all NA vals and we will skip these after marking NA
        #Don't need to add an index in the plot list
        #might be smarter to just not add these in the future but don't want to deal with the indexing now
        if(all(is.na(y))){
            gr.df[j, "iso_name"] = iso_name
            gr.df[j, "gr.est"] = NA
            gr.df[j, "gr.p"] = NA
            gr.df[j, "gr.rsq"] = NA
            gr.df[j, "n_obs"] = 0
            next
        }
        
        #run reg
        temp_lm = lm(y ~ x)
        x_est = summary(temp_lm)[[4]][2,1]
        x_p = summary(temp_lm)[[4]][2,4]
        x_rsq = summary(temp_lm)[[8]]
        #add stats to table
        gr.df[j, "iso_name"] = iso_name
        gr.df[j, "gr.est"] = x_est
        gr.df[j, "gr.p"] = x_p
        gr.df[j, "gr.rsq"] = x_rsq
        gr.df[j, "n_obs"] = sum( !is.na(y))
        
        
        #plot
        p_title = paste(iso_name, "set", iso_set, "rep", iso_rep, "temp", iso_temp, sep = " ")
        p = ggplot(
            data.frame(days = x, mm = y),
            aes(days, mm)
        ) +
            geom_point() +
            geom_smooth(method = "lm") +
            theme_bw() +
            labs(title = p_title)
        
        #cat setRepTemp to index list
        setRepTemp = paste(
            set_rep_temp.list[[i]][1],
            set_rep_temp.list[[i]][2],
            set_rep_temp.list[[i]][3],
            sep = "."
        )
        
        plots.list[[iso_name]][[setRepTemp]] = p
    }
    gr_stat.list[[i]] = gr.df
}
gr_stat.df = bind_rows(gr_stat.list)
nrow(gr_stat.df)
write.table(gr_stat.df, file = "data/summary_tables/growth_rate_stats.txt", col.names = T, row.names = F, sep = "\t", quote = F)

############################
#Plots and data checks

gr_stat.df %>% filter(gr.rsq < 0.97)
gr_stat.df %>% filter(n_obs < 5)

#there was a warning about perfect fit. Check for R^2 of 1 or close to
gr_stat.df %>% filter(gr.rsq > 0.9999)
#There are 9 isolates but none have perfect fit

#n_obs v rsq
p = ggplot(gr_stat.df, aes(x = n_obs, y = gr.rsq)) +
    geom_point() +
    geom_hline(yintercept = 0.97) +
    theme_bw()
p

pdf("figures/n_obs_v_rsq.pdf")
p
dev.off()

#plot 95 <= growth curves < 97% r-squared
#plot 90 <= growth curves < 95% r-squared
# plot growth curves < 90% r-squared
lt_97 = gr_stat.df %>% filter(gr.rsq < 0.97 & gr.rsq >= 0.95)
lt_95 = gr_stat.df %>% filter(gr.rsq < 0.95 & gr.rsq >= 0.90)
lt_90 = gr_stat.df %>% filter(gr.rsq < 0.90 )

nrow(lt_97)
nrow(lt_95)
nrow(lt_90)

lt_97$temp %>% unique()
lt_95$temp %>% unique()
lt_90$temp %>% unique()

pdf("figures/lt_97_rsq.pdf")
for(i in 1:nrow(lt_97)){
    setRepTemp = paste(lt_97[i, "set"], lt_97[i, "rep"], lt_97[i, "temp"], sep = ".")
    iso_name = lt_97[i, "iso_name"]
    print( plots.list[[iso_name]][[setRepTemp]] )
}
dev.off()

pdf("figures/lt_95_rsq.pdf")
for(i in 1:nrow(lt_95)){
    setRepTemp = paste(lt_95[i, "set"], lt_95[i, "rep"], lt_95[i, "temp"], sep = ".")
    iso_name = lt_95[i, "iso_name"]
    print( plots.list[[iso_name]][[setRepTemp]] )
}
dev.off()

pdf("figures/lt_90_rsq.pdf")
for(i in 1:nrow(lt_90)){
    setRepTemp = paste(lt_90[i, "set"], lt_90[i, "rep"], lt_90[i, "temp"], sep = ".")
    iso_name = lt_90[i, "iso_name"]
    print( plots.list[[iso_name]][[setRepTemp]] )
}
dev.off()


#########
#Check number reps per isolate x temp
n_reps = gr_stat.df %>%
    filter(!is.na(gr.est)) %>%
    group_by(iso_name, temp, set) %>%
    summarize(n_reps = n())
head(n_reps)

filter(n_reps, n_reps < 3 & set == 1) %>% print(n = Inf)
#Only 5 and 35 C (2 reps currently run)
filter(n_reps, n_reps < 3 & set == 2) %>% print(n = Inf)
#Only 5 and 35 C (2 reps currently run)
filter(n_reps, n_reps < 3) %>% print(n = Inf)

# write.csv(filter(n_reps, n_reps < 3) %>% print(n = Inf), "data/blank_data_sheets/redo_list_05202024.csv", row.names = F)

#####
#print all plots to giant pdf
pdf("figures/gr_lms.pdf")

#lapply(plots.list, print) #this prints a ton of extra output to the console

for(i in 1:length(plots.list)){
    for(j in 1:length(plots.list[[i]])){
        print(plots.list[[i]][[j]])
    }
}

dev.off()
