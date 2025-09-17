library(dplyr)
library(data.table)

#############################
#FUNS
#############################

#function for nearest neighbor
#we will have already searched for SNPs *within* genes so the we can then search for start sites > then pos else stop sites < pos for the nearest neighbor match
#genes down stream of SNPs will have positive distance and upstream will have negative distance
#i.e., start - SNPpos OR stop - SNPpos
nearest_gene = function(snpPos = NULL, gff = NULL){
    #new col for distances
    gff$dist = vector(mode = "numeric", length = nrow(gff))
    #calculate distances
    gff$dist = ifelse(gff$start > snpPos, gff$start - snpPos, gff$stop - snpPos)
    #filter to minimum absolute distance and return the resulting df
    #dplyr way
    #return( filter(gff, abs(dist) == min(abs(gff$dist))) )
    #DT way
    return(
        gff[
            abs(dist) == min(abs(gff$dist))
        ]
    )
}

############################################
#nearest neighbor test routine START

#foo = data.frame(start = seq(1:10), stop = 1:10, dist = vector(mode = "numeric", length = 10))

#foo
#snpPos = 5.5
#foo$dist = ifelse(foo$start > snpPos, foo$start - snpPos, foo$stop - snpPos)
#foo

#filter(foo, abs(dist) == min(abs(foo$dist)))

#foo2 = data.frame(start = seq(1:10), stop = 1:10)
#foo3 = nearest_gene(snpPos, foo2)
#foo3

#test routine END
############################################


##########################
#Import data
##########################

#import gff (gene positions) contig name, gene ID, start site, stop site
gff = read.csv("data/genotype_dat/Neonectria_faginata.funannotate.mRNA.gff3", sep = "\t", header = F)
head(gff)
gff = gff[,c(1,9,4,5)]
colnames(gff) = c("contig", "geneID", "start", "stop")
head(gff)

#import LFMM results
lfmm_results = read.table("data/LFMM/grpc1.lfmm.txt", header = T)
colnames(lfmm_results)
#filter for SNPs with sig relationship to variable
############
#USE MAN ADJUSTED P 
lfmm_results.sig = lfmm_results %>% filter(FDR.sig.man == "sig") %>% select(scaffold, position)
nrow(lfmm_results.sig)
head(lfmm_results.sig)
lfmm_results.sig$scaffold = sub("_pilon", "", lfmm_results.sig$scaffold)

#################################
#New cols for reporting matches
#################################

#new col in gff table for snp counts
gff$sig_SNP_count = rep(0, nrow(gff))
#new cols in SNP table for whether a gene was found for a snp and type of relationship (no.match, in.gene, nearest), gene id of match, and distance of match (if nearest neighbor)
lfmm_results.sig$match.type = vector(mode = "character", length = nrow(lfmm_results.sig))
lfmm_results.sig$geneID = vector(mode = "character", length = nrow(lfmm_results.sig))
lfmm_results.sig$distance = vector(mode = "numeric", length = nrow(lfmm_results.sig))
lfmm_results.sig$indel = vector(mode = "character", length = nrow(lfmm_results.sig))

#set to data.table to using rolling count update
gff.dt = as.data.table(gff)

#################################################

#Loop through SNPs
#Need to add gene name matches to SNP table (also deal with case of >1 match)
#Also add matched SNP count to genes 
for(i in 1:nrow(lfmm_results.sig)){
    
    #filter for a gene that contains the SNP within start and stop positions
    temp = gff.dt[
        contig == lfmm_results.sig$scaffold[i] &
        start <=  lfmm_results.sig$position[i] &
        stop >= lfmm_results.sig$position[i]
    ]
    #check if there is a matching gene. If not go to next
    if(nrow(temp) == 1){ #need the if in cases where there is no in-gene match
        gff.dt[ #the first lines are dt conditionals
            contig == lfmm_results.sig$scaffold[i] &
                start <=  lfmm_results.sig$position[i] &
                stop >= lfmm_results.sig$position[i],
            #rolling count update
            sig_SNP_count := sig_SNP_count + 1
        ] 
    
        lfmm_results.sig$match.type[i] = "in.gene"
        lfmm_results.sig$geneID[i] = temp$geneID
        lfmm_results.sig$distance[i] = 0
    }else{ #if there is no within gene match filter to the nearest gene(s)
        #filter to the correct scaffold
        temp2 = gff.dt[
            contig == lfmm_results.sig$scaffold[i] 
        ] 
        #calculate distance to genes and return the nearest neighbor(s)
        temp2.match = nearest_gene(lfmm_results.sig$position[i], temp2)
        
        #Need to account for possible ties in distance. First look for single match
        if(nrow(temp2.match) == 1){
            print(nrow(temp2.match))
            gff.dt[
                geneID == temp2.match$geneID, #conditional filter to gene ID
                sig_SNP_count := sig_SNP_count + 1 #update count
            ]  
      
            lfmm_results.sig$match.type[i] = "nearest"
            lfmm_results.sig$geneID[i] = temp2.match$geneID
            lfmm_results.sig$distance[i] = temp2.match$dist
        }else if(nrow(temp2.match) > 0){ #if there are no matches we will update below
            print(nrow(temp2.match))
            gff.dt[
                geneID == temp2.match$geneID, #conditional filter to gene ID
                sig_SNP_count := sig_SNP_count + 1 #update count
            ]  
            
            #lfmm_results.sig$match.type[i] = "nearest.multiple"
            #lfmm_results.sig$geneID = temp2.match$geneID
            #lfmm_results.sig$distance = temp2.match$dist
        }else{
            lfmm_results.sig$match.type[i] = "no.match"
            lfmm_results.sig$geneID[i] = "no.match"
            lfmm_results.sig$distance[i] = 0
        }
    }#end alternative to within gene
}#end for loop
            
lfmm_results.sig$match.type %>% unique
lfmm_results.sig$distance %>% abs %>% mean #2032.842
lfmm_results.sig$distance %>% abs %>% median #799
lfmm_results.sig %>% filter(abs(distance) > 0) %>% pull(distance) %>% abs %>% range #109 9178
lfmm_results.sig %>% filter(match.type == "in.gene") %>% nrow #6
lfmm_results.sig %>% nrow #6 of 19 in gene

gff.dt %>% filter(sig_SNP_count > 0) #18 genes have SNPs
gff.dt %>% filter(sig_SNP_count > 1) #1 have at least 2 SNPs
nrow(gff.dt)

#######################
#plots of snp distance and SNP accumulation
#######################
library(ggplot2)
p1 = ggplot(lfmm_results.sig, aes(x = distance)) +
    geom_histogram(breaks = seq(-15000, 15000, 250)) +
    #geom_histogram(breaks = c(-15000,-12500,-10000,-7500,-5000, -2500, -1000, -500,-250, -5, 5, 250, 500, 1000, 2500, 5000, 7500, 10000, 1250, 1500)) +
    theme_bw() +
    labs(x = "Gene distance from SNP (bp)", y = "Number SNPs (bin width 250 bp)")
p1
pdf("figures/LFMM/SNP_distance.grpc1.pdf")
p1
dev.off()

#x = reorder(geneID, -sig_SNP_count), 

p2 = ggplot(
        gff.dt %>% filter(sig_SNP_count > 0), 
        aes(x = sig_SNP_count)
    ) +
    geom_bar() +
    labs(x = "SNPs per gene", y = "No. genes") +
    theme_bw()

p2
pdf("figures/LFMM/SNP_per_gene.grpc1.pdf")
p2
dev.off()

ggplot(
    gff.dt %>% filter(sig_SNP_count > 0), 
    aes(x = reorder(geneID, -sig_SNP_count), y = sig_SNP_count)
) +
    geom_point() +
    labs(x = "Gene", y = "No. SNPs") +
    theme_bw() +
    theme(
        axis.text.x = element_blank()
    )

##############################
#Write results tables
##############################
write.table(gff.dt %>% filter(sig_SNP_count > 0), "data/LFMM/annotations/grpc1.gene_SNP_hits.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(lfmm_results.sig, "data/LFMM/annotations/grpc1.SNPs.gene_found.nearest_neighbors.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(gff.dt %>% filter(sig_SNP_count > 0) %>% select(geneID), "data/LFMM/annotations/grpc1.geneIDs.nearest_neighbors.txt", col.names = F, row.names = F, sep = "\t", quote = F)
