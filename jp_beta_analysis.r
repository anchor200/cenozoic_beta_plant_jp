library(ClassDiscovery)
library(maps)
library(ggplot2)
library(reshape2)
library(betapart)
library(vegan)
library(geosphere)

# reading geographic information
regions <- read.delim("data/region_mask.tsv", sep='\t', header=F, stringsAsFactor=FALSE)
colnames(regions) <- 120:149
rownames(regions) <- 45:20

seas <- read.delim("data/sep_sea.tsv", sep='\t', header=F, stringsAsFactor=FALSE)
colnames(seas) <- 120:149
rownames(seas) <- 45:20

islands <- read.delim("data/sep_island.tsv", sep='\t', header=F, stringsAsFactor=FALSE)
colnames(islands) <- 120:149
rownames(islands) <- 45:20

culdata <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
colnames(culdata) <- c("lat", "lon", "region", "sea", "island")

for (i in 120:149){
    for (j in 20:45){
        temp <- regions[as.character(j), as.character(i)]
        temps <- seas[as.character(j), as.character(i)]
        tempi <- islands[as.character(j), as.character(i)]
        if (temp!=-1){
            tempRow <- data.frame("lat"=j, "lon"=i, "region"=temp, "sea"=temps, "island"=tempi)
            culdata <- rbind(culdata, tempRow)
        }
    }
}

# naming regions
culdata$Okhotsk <- 0
culdata$Okhotsk[culdata$sea==0] <- 1
culdata$Pacific1 <- 0
culdata$Pacific1[culdata$sea==2] <- 1
culdata$Pacific2 <- 0
culdata$Pacific2[culdata$sea==3] <- 1
culdata$Pacific3 <- 0
culdata$Pacific2[culdata$sea==6] <- 1
culdata$Japan <- 0
culdata$Japan[culdata$sea==4] <- 1
culdata$China <- 0
culdata$China[culdata$sea==5] <- 1

culdata$Hokkaido <- 0
culdata$Hokkaido[culdata$island==0] <- 1
culdata$Honshu <- 0
culdata$Honshu[culdata$island==1] <- 1
culdata$Shikoku <- 0
culdata$Shikoku[culdata$island==3] <- 1
culdata$Kushu <- 0
culdata$Kushu[culdata$island==2] <- 1
culdata$Okinawa <- 0
culdata$Okinawa[culdata$island==4] <- 1
culdata$Okinawa2 <- 0
culdata$Okinawa2[culdata$island==5] <- 1

culdata$key <- paste(culdata$lat, culdata$lon, sep="_")
culdata <- culdata[, colnames(culdata)!="sea" & colnames(culdata)!="lat" & colnames(culdata)!="lon" & colnames(culdata)!="island"]
culdata_resh$region_and_key <- paste(culdata_resh$key, culdata_resh$region, sep="_")


# devision of Hokkaido
culdata_resh_Hokkaido <- subset(culdata_resh, culdata_resh$Hokkaido==1)
culdata_resh_kai <- culdata_resh_Hokkaido[, colnames(culdata_resh_Hokkaido)!="key" & colnames(culdata_resh_Hokkaido)!="region" & colnames(culdata_resh_Hokkaido)!="region_and_key"]
culdata_resh_mat <- data.frame(t(as.matrix(culdata_resh_kai)))
colnames(culdata_resh_mat) <- culdata_resh_Hokkaido$key

site_dist <- distanceMatrix(culdata_resh_mat, metric="euclidian")
site.cluster <- hclust(site_dist)
d <- as.dendrogram(site.cluster)
site.dend <- as.dendrogram(d)
colnames(culdata_resh_mat) <- culdata_resh_Hokkaido$region_and_key

clusters<-cutree(site.cluster,k=3)
## clusters<-cutree(site.cluster,k=2) ## for supplementary analysis

cluster_resh <- data.frame(clusters)
rownames(cluster_resh) <- culdata_resh_Hokkaido$region_and_key
cluster_resh$lat <- 0
cluster_resh$lon <- 0
for (i in 1:nrow(cluster_resh)){
    cluster_resh$lat[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][1])
    cluster_resh$lon[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][2])
}
cluster_resh_n <- cluster_resh

# division of Honshu, Shikoku, Kyushu
culdata_resh_Honshu <- subset(culdata_resh, culdata_resh$Hokkaido==0 & culdata_resh$Okinawa==0& culdata_resh$Okinawa2==0)
culdata_resh_Honshu <- subset(culdata_resh_Honshu, !(culdata_resh_Honshu$lat==41 & lon==139))
culdata_resh_kai <- culdata_resh_Honshu[, colnames(culdata_resh_Honshu)!="key" & colnames(culdata_resh_Honshu)!="region" & colnames(culdata_resh_Honshu)!="region_and_key"]
culdata_resh_mat <- data.frame(t(as.matrix(culdata_resh_kai)))
colnames(culdata_resh_mat) <- culdata_resh_Honshu$key

site_dist <- distanceMatrix(culdata_resh_mat, metric="euclidian")
site.cluster <- hclust(site_dist)
d <- as.dendrogram(site.cluster)
site.dend <- as.dendrogram(d)
colnames(culdata_resh_mat) <- culdata_resh_Honshu$region_and_key

clusters<-cutree(site.cluster,k=14)
## clusters<-cutree(site.cluster,k=10) ## for supplementary analysis

cluster_resh <- data.frame(clusters)
rownames(cluster_resh) <- culdata_resh_Honshu$region_and_key
cluster_resh$clusters <- cluster_resh$clusters + 3
cluster_resh$lat <- 0
cluster_resh$lon <- 0
for (i in 1:nrow(cluster_resh)){
    cluster_resh$lat[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][1])
    cluster_resh$lon[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][2])
}
cluster_resh_c <- cluster_resh

# merge & write division table
site_division <- rbind(cluster_resh_n, cluster_resh_c)
write.table(site_division, "data/out/site_division.tsv", sep="\t", quote=F, row.names=F)
## write.table(site_division, "data/out/site_division_suppl.tsv", sep="\t", quote=F, row.names=F) ## for supplementary analysis

geo.ages <- c("Oligocene", "Miocene", "Pliocene", "Pleistocene", "LastGlacial", "Holocene", "Present")

# Tertiary genus occurrence
data_resh<-read.delim("data/Japan_Tertiary.csv", sep=',', header=T, stringsAsFactor=FALSE)
data_resh$Genus <- tolower(data_resh$Genus)
data_resh$lat_floor <- floor(data_resh$Latitude)
data_resh$long_floor <- floor(data_resh$Longitude)
data_resh$lat_long <- paste(data_resh$lat_floor, data_resh$long_floor, sep="_")

# Quartenary genus occurrence
data_qua<-read.delim("data/Japan_Quaternary.csv", sep=',', header=T, stringsAsFactor=FALSE)
data_qua$Genus <- tolower(data_qua$genus)
data_qua$Latitude <- data_qua$lat
data_qua$Longitude <- data_qua$lon
data_qua$lat_floor <- floor(data_qua$Latitude)
data_qua$long_floor <- floor(data_qua$Longitude)
data_qua$lat_long <- paste(data_qua$lat_floor, data_qua$long_floor, sep="_")

# Current genus occurrence
dataresent_2d<-readRDS("data/present_tree_mesh2.rds")
dataresent_2d$Latitude <- dataresent_2d$lat
dataresent_2d$Longitude <- dataresent_2d$lon
dataresent_2d$lat_floor <- floor(dataresent_2d$Latitude)
dataresent_2d$long_floor <- floor(dataresent_2d$Longitude)
dataresent_2d$Genus <- tolower(dataresent_2d$genus)
dataresent_2d$lat_long <- paste(dataresent_2d$lat_floor, dataresent_2d$long_floor, sep="_")


# Data description
print(c("Tertiary data points", nrow(data_resh)))
print(c("Tertiary sites", nrow(table(data_resh$Plot.name))))
print(c("Tertiary n genus", nrow(table(data_resh$Genus))))

print(c("Quartenary data points", nrow(data_qua)))
print(c("Quartenary sites", nrow(table(data_qua$chimei))))
print(c("Quartenary n genus", nrow(table(data_qua$Genus))))

table(data_qua$花粉vs大型化石)

# Genus list
genus_list <-union(union(union(rownames(table(dataresent_2d$Genus)), rownames(table(data_qua$Genus))), rownames(table(data_resh$Genus))), rownames(table(dataresent_2d$Genus)))
nGen <- length(genus_list)
genus_list_ter <-rownames(table(data_resh$Genus))
nGen_ter <- length(genus_list_ter)


# plot of fossile sites

## options(repr.plot.width=4, repr.plot.height=4)
plot_co_places <- function(cene){

    data_cene <- subset(data_resh, data_resh$Geological.period==cene)    
    
    places <- rownames(table(data_cene$Plot.name))
    n_places <- nrow(table(data_cene$Plot.name))

    world.map <- map_data ("world")
    japan <- world.map[world.map$long >= 120 & world.map$long < 150 & world.map$lat >= 20 & world.map$lat < 46 & world.map$region=="Japan",]
    g <- ggplot(japan, aes(x = long, y = lat)) + geom_path(aes(group = group), color="blue")+
      layer(
        data=data_cene, 
        mapping=aes(x=Longitude, y=Latitude), 
        geom="point", 
        stat="identity", 
        position="identity"
      ) + ggtitle(paste(cene, "(",as.character(n_places), "sites)"))+theme_bw()+ theme(text = element_text(size = 24))
    return(g)
}

plot_co_places_qua <- function(cene){

    data_cene <- subset(data_qua, data_qua$Geological.period==cene)    
    
    places <- rownames(table(data_cene$Latitude))
    n_places <- nrow(table(data_cene$Latitude))
    world.map <- map_data ("world")
    japan <- world.map[world.map$long >= 120 & world.map$long < 150 & world.map$lat >= 20 & world.map$lat < 46 & world.map$region=="Japan",]
    g <- ggplot(japan, aes(x = long, y = lat)) + geom_path(aes(group = group), color="blue")+
      layer(
        data=data_cene, 
        mapping=aes(x=Longitude, y=Latitude), 
        geom="point", 
        stat="identity", 
        position="identity"
      ) + ggtitle(paste(cene, "(", as.character(n_places), "sites)"))+theme_bw()+ theme(text = element_text(size = 24))
    return(g)

}

options(repr.plot.width=24, repr.plot.height=16)
p1 <- plot_co_places("Oligocene")
p2 <- plot_co_places("Miocene")
p3 <- plot_co_places("Pliocene")
p4 <- plot_co_places_qua("Pleistocene")
p5 <- plot_co_places_qua("Last glacial period")
p6 <- plot_co_places_qua("Holocene")

gridExtra::grid.arrange( p1, p2, p3, p4, p5, p6, nrow = 2)

# read a sites data
sepfile <- "data/out/site_division_matrix.tsv"
regions <- read.delim(sepfile, sep='\t', header=F, stringsAsFactor=FALSE)

colnames(regions) <- 120:149
rownames(regions) <- 45:20
regions$lat <- rownames(regions)
regions_resh <- melt(regions, id.vars=c("lat"), variable.name="long")
regions_resh$lat_long <- paste(regions_resh$lat, regions_resh$long, sep="_")
regions_resh$lat <- as.numeric(regions_resh$lat)
regions_resh$long <- as.numeric(as.character(regions_resh$long))

## lat / lon values represent lower boundary of grid.


n_regions = 100 # just a big number for produce a loop
regions = 0:(n_regions-1)

regions_info <- data.frame(matrix(0, nrow=n_regions))
colnames <- c("region_id")
regions_info$region_id <- regions
regions_info <- regions_info[-1]
##  centering grid axes
for (id in regions){
    regions_info$lat_cent[regions_info$region_id==id] <- mean(subset(regions_resh, regions_resh$value==id)$lat) + 0.5
    regions_info$long_cent[regions_info$region_id==id] <- mean(subset(regions_resh, regions_resh$value==id)$long) + 0.5
}

rownames(regions_info) <- regions_info$region_id
regions_info <- regions_info[-1]

# extended default functions
merge2 <- function(dfs, ...)
{
　 base <- dfs[1]
　 lapply(dfs[-1], function(i) base <<- merge(base, i, ...)) # [1]
  return(base)
}
mean_ = function(x){
    return(mean(x, na.rm=TRUE))
}
sd_ = function(x){
    return(sd(x, na.rm=TRUE))
}

colname_correction_region <- function(cene){

    pa_by_region <- get(paste(cene, "_pa_by_region", sep=""))
    names <- colnames(pa_by_region)
    for (n in 1:length(names)){
        names[n] <- paste(cene, names[n], sep="")
    }
    colnames(pa_by_region) <- names
    temp <- as.data.frame(pa_by_region)
    temp$genus <- rownames(temp)
    return(temp)
}

whole_samples_region <- merge2(list(
    colname_correction_region("Oligocene"),
    colname_correction_region("Miocene"),
    colname_correction_region("Pliocene"),
    colname_correction_region("Pleistocene"),
    colname_correction_region("LastGlacial"),
    colname_correction_region("Holocene"),
    colname_correction_region("Present")
),by=c("genus"))

## data reshape for analyses
rownames(whole_samples_region) <- whole_samples_region$genus
for_ana_region <- whole_samples_region[-1]

## sample completeness calculation
pool_1 <- t(for_ana_region)
pool_1 <- pool_1[apply(pool_1, 1, sum)>1,]
pool_ <- gsub("\\d", "", rownames(pool_1))
pool <- specpool(pool_1, pool_)
pool$coverage <- pool$Species / pool$chao

ppp <- pool_1[1:4,]


## filtering with n of genus par subregion
CUT <- 10
for_ana_region[for_ana_region > 1] <- 1
for_ana_region <- for_ana_region[apply(for_ana_region, 2, sum)>CUT]


#Tertiary data processing func
plot_co_region_matrix <- function(cene){

    data_cene <- subset(data_resh, data_resh$Geological.period==cene) 
    
    
    temp_aggr <- matrix(0, nrow=nGen, ncol=n_regions)
    rownames(temp_aggr) <- genus_list
    colnames(temp_aggr) <- regions

    for (reg in regions){
        place = subset(regions_resh, regions_resh$value==reg)$lat_long
        for (p in place){
            temp <- subset(data_cene, data_cene$lat_long==p)
            for (sp in temp$Genus){
                temp_aggr[sp, as.character(reg)] = temp_aggr[sp, as.character(reg)]  + 1
            }
        }
    }
    return(temp_aggr)

}

#Quartenary data processing func
plot_co_region_matrix_qua <- function(cene){

    data_cene <- subset(data_qua, data_qua$Geological.period==cene)    
    
    temp_aggr <- matrix(0, nrow=nGen, ncol=n_regions)
    rownames(temp_aggr) <- genus_list
    colnames(temp_aggr) <- regions

    for (reg in regions){
        place = subset(regions_resh, regions_resh$value==reg)$lat_long
        for (p in place){
            temp <- subset(data_cene, data_cene$lat_long==p)
            for (sp in temp$Genus){
                temp_aggr[sp, as.character(reg)] = temp_aggr[sp, as.character(reg)]  + 1
            }
        }
    }
    return(temp_aggr)

}

Oligocene_pa_by_region <- plot_co_region_matrix("Oligocene")
Miocene_pa_by_region <- plot_co_region_matrix("Miocene")
Pliocene_pa_by_region <- plot_co_region_matrix("Pliocene")

Pleistocene_pa_by_region <- plot_co_region_matrix_qua("Pleistocene")
LastGlacial_pa_by_region <- plot_co_region_matrix_qua("Last glacial period")
Holocene_pa_by_region <- plot_co_region_matrix_qua("Holocene")

#Current data processing
    data_cene <-   dataresent_2d
    temp_aggr <- matrix(0, nrow=nGen, ncol=n_regions)
    rownames(temp_aggr) <- genus_list
    colnames(temp_aggr) <- regions
    for (reg in regions){
        place = subset(regions_resh, regions_resh$value==reg)$lat_long
        for (p in place){
            temp <- subset(data_cene, data_cene$lat_long==p)
            for (sp in temp$Genus){
                temp_aggr[sp, as.character(reg)] = temp_aggr[sp, as.character(reg)]  +1
            }
        }
    }


Present_pa_by_region <- temp_aggr


# function for beta div calc

beta_order <- function(res, info=0){
    size <- ncol(as.matrix(res$beta.sim))
    n_pairs <- size * (size - 1) / 2
    output <- matrix(0, ncol=6, nrow=n_pairs)
    e_regions <- union(rownames(as.matrix(res$beta.sim)), colnames(as.matrix(res$beta.sim)))
    pairs <- combn(e_regions, m=2)        
    colnames(output) <- c("s_lat", "n_lat", "beta.sim", "beta.sne", "beta.sor", "info")

    for (i in 1:(n_pairs)){
        output[i,"s_lat"] <- gsub("\\D", "", pairs[,i][1])
        output[i,"n_lat"] <- gsub("\\D", "", pairs[,i][2])
        output[i,"beta.sim"] <- as.matrix(res$beta.sim)[pairs[,i][1],pairs[,i][2]]
        output[i,"beta.sne"] <- as.matrix(res$beta.sne)[pairs[,i][1],pairs[,i][2]]
        output[i,"beta.sor"] <- as.matrix(res$beta.sor)[pairs[,i][1],pairs[,i][2]]
        output[i,"info"] <- info
    }
    output <- as.data.frame(output)
    output$beta.sim <- as.numeric(output$beta.sim)
    output$beta.sne <- as.numeric(output$beta.sne)
    output$beta.sor <- as.numeric(output$beta.sor)
    return(output)
}

# function for regression
beta_show4 <- function(lat, ana=for_ana){
    res <- beta.pair(t(as.matrix(select(ana, contains(lat)))))
    beta_order(res, lat)
}

beta.regr <- function(cene, data_aggr, method, clim=FALSE){
    if(clim==FALSE){
    da <- subset(data_aggr, data_aggr$info == cene)
    da$temp <- da$beta.sim
    
    if (method == "pow"){
        da$sp_dist <- log(da$sp_dist)
        res <- glm(temp+.Machine$double.eps ~ sp_dist , family = gaussian(log), data=da)
    }else{
        res <- glm(temp+.Machine$double.eps ~ sp_dist , family = gaussian(log), data=da)
    }
    res.temp <- as.data.frame(coef(summary(res)))
    res.temp$geo_peri <- cene
    
    return(list("coef"=res.temp, "pseudo_r2"=1 - res$deviance / res$null.deviance, "AIC"=AIC(res), "beta.mean"=mean_(da$beta.sim), "beta.sd"=sd_(da$beta.sim), "result"=res))

    }
    else{

    da <- subset(data_aggr, data_aggr$info == cene)
    if (cene == "LastGlacial"){
        da$clim_dist <- da[,"clim_dist_LGP"]
    }
    else if(cene == "Present"){
        da$clim_dist <- da[,"clim_dist_Present"]
    }
    else{
        da$clim_dist <- da[,paste("clim_dist_", gsub("cene", "", cene), sep="")]        
    }


    if (method == "pow"){
        da$clim_dist <- log(da$clim_dist + .Machine$double.eps)
        res <- glm(1 - beta.sim+.Machine$double.eps ~ clim_dist , family = gaussian(log), data=da)        
    }else{
            res <- glm(1 - beta.sim+.Machine$double.eps ~ clim_dist , family = gaussian(log), data=da)
    }
    res.temp <- as.data.frame(coef(summary(res)))
    res.temp$geo_peri <- cene
    
    return(list("coef"=res.temp, "pseudo_r2"=1 - res$deviance / res$null.deviance, "AIC"=AIC(res), "beta.mean"=mean_(da$beta.sim), "beta.sd"=sd_(da$beta.sim), "result"=res))

    }

}



# regression result reshape
beta.calc.by.ages <- function(method, xseq, ret, clim=FALSE){
    if (xseq == "clim"){
        return
    }
    else{
        dd <- plot_data_region_dist
    }

    
    if (ret == "AIC"){
        temp.ret <- rbind(
        beta.regr("Oligocene", dd, method, clim)$AIC,
        beta.regr("Miocene", dd, method, clim)$AIC,
        beta.regr("Pliocene", dd, method, clim)$AIC,
        beta.regr("Pleistocene", dd, method, clim)$AIC,
        beta.regr("LastGlacial", dd, method, clim)$AIC,
        beta.regr("Holocene", dd, method, clim)$AIC,
        beta.regr("Present", dd, method, clim)$AIC
     )
    }else if (ret == "pseudo_r2"){
        temp.ret <- rbind(
        beta.regr("Oligocene", dd, method, clim)$pseudo_r2,
        beta.regr("Miocene", dd, method, clim)$pseudo_r2,
        beta.regr("Pliocene", dd, method, clim)$pseudo_r2,
        beta.regr("Pleistocene", dd, method, clim)$pseudo_r2,
        beta.regr("LastGlacial", dd, method, clim)$pseudo_r2,
        beta.regr("Holocene", dd, method, clim)$pseudo_r2,
        beta.regr("Present", dd, method, clim)$pseudo_r2
         )
    }else{
        
        temp.ret <- rbind(
        beta.regr("Oligocene", dd, method, clim)$coef,
        beta.regr("Miocene", dd, method, clim)$coef,
        beta.regr("Pliocene", dd, method, clim)$coef,
        beta.regr("Pleistocene", dd, method, clim)$coef,
        beta.regr("LastGlacial", dd, method, clim)$coef,
        beta.regr("Holocene", dd, method, clim)$coef,
        beta.regr("Present", dd, method, clim)$coef
     )
    }

return(temp.ret)

}



# beta diversity calculation
plot_data_region_dist <- rbind(
beta_show4("Oligocene", for_ana_region),
beta_show4("Miocene", for_ana_region),
beta_show4("Pliocene", for_ana_region),
beta_show4("Pleistocene", for_ana_region),
beta_show4("LastGlacial", for_ana_region),
beta_show4("Holocene", for_ana_region),
beta_show4("Present", for_ana_region)
) 

plot_data_region_dist$info <- factor(plot_data_region_dist$info, levels=geo.ages)
colnames(plot_data_region_dist) <- c("regionA", "regionB", "beta.sim", "beta.sne",  "beta.sor", "info")
plot_data_region_dist$regions <- paste(plot_data_region_dist$regionA, plot_data_region_dist$regionB, sep="_")
plot_data_region_dist$regionA_cent_lat<- regions_info[plot_data_region_dist$regionA,1]
plot_data_region_dist$regionA_cent_long<- regions_info[plot_data_region_dist$regionA,2]
plot_data_region_dist$regionB_cent_lat<- regions_info[plot_data_region_dist$regionB,1]
plot_data_region_dist$regionB_cent_long<- regions_info[plot_data_region_dist$regionB,2]
plot_data_region_dist$dist <- 0

for (k in 1:nrow(plot_data_region_dist)){
    plot_data_region_dist$dist[k] <- distGeo(c(plot_data_region_dist$regionA_cent_long[k], plot_data_region_dist$regionA_cent_lat[k]), 
                                             c(plot_data_region_dist$regionB_cent_long[k], plot_data_region_dist$regionB_cent_lat[k])) / 1000
}

plot_data_region_dist$label <- paste(plot_data_region_dist$regionA, plot_data_region_dist$regionB, sep="_")
plot_data_region_dist$info <- factor(plot_data_region_dist$info, levels=geo.ages)
plot_data_region_dist$sp_dist <- plot_data_region_dist$dist


# spatial distance-decay
## negative exp
neg_exp_coefs_distonly_region <- beta.calc.by.ages("neg", "space", "ret")
neg_exp_coefs_distonly_region$predictor <- gsub("\\d", "", rownames(neg_exp_coefs_distonly_region))
colnames(neg_exp_coefs_distonly_region) <- gsub("\\s", "",colnames(neg_exp_coefs_distonly_region))
neg_exp_coefs_distonly_region

## power law
pow_coefs_distonly_region <- beta.calc.by.ages("pow", "space", "ret")

pow_coefs_distonly_region$predictor <- gsub("\\d", "", rownames(pow_coefs_distonly_region))
colnames(pow_coefs_distonly_region) <- gsub("\\s", "",colnames(pow_coefs_distonly_region))


# fit measures
neg_exp_aic <- beta.calc.by.ages("neg", "space", "AIC", clim=FALSE)
neg_exp_r2 <- beta.calc.by.ages("neg", "space", "pseudo_r2", clim=FALSE)
neg_exp_fit <- data.frame("info"=geo.ages, "pseudo_r2"=neg_exp_r2, "AIC"=neg_exp_aic)
neg_exp_fit

# climatic distance-decay
paleoclim<-read.delim("data/paleoclim_degree.csv", sep=',', header=T, stringsAsFactor=FALSE)
paleoclim$int_gridkey <- paste(floor(paleoclim$y), floor(paleoclim$x), sep="_")
paleoclim_resh <- paleoclim[c("int_gridkey", "Oligo", "Mio", "Plio", "Pleist", "LGP", "Holo", "current")]
temp_merged <- merge(regions_resh,paleoclim_resh,by.x="lat_long", by.y="int_gridkey", all=T)
clim_region <- subset(temp_merged, temp_merged$value>=0)

mean__ <- function(x){
    return(list("mean"=mean_(x), "sd"=sd_(x)))
}
clim_resh <- aggregate(list(
               "Oligo"=clim_region$Oligo,
               "Mio"=clim_region$Mio,
               "Plio"=clim_region$Plio,
               "Pleist"=clim_region$Pleist,
               "LGP"=clim_region$LGP,
               "Holo"=clim_region$Holo,
               "Present"=clim_region$current),by=list("region"=clim_region$value), mean__)

# compile climate data
for (k in 1:nrow(plot_data_region_dist)){
    plot_data_region_dist$clim_dist_Oligo[k] <- abs(clim_resh$Oligo[[as.numeric(plot_data_region_dist$regionA[k])+1,1]]-
    clim_resh$Oligo[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Mio[k] <- abs(clim_resh$Mio[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Mio[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Plio[k] <- abs(clim_resh$Plio[[as.numeric(plot_data_region_dist$regionA[k])+1,1]]-
    clim_resh$Plio[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Pleisto[k] <- abs(clim_resh$Pleist[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Pleist[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_LGP[k] <- abs(clim_resh$LGP[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$LGP[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Holo[k] <- abs(clim_resh$Holo[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Holo[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Present[k] <- abs(clim_resh$Present[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Present[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    
}

## negative exp
neg_exp_coefs_clim_region <- beta.calc.by.ages("neg", "space", "ret", clim=TRUE)
neg_exp_coefs_clim_region$predictor <- gsub("\\d", "", rownames(neg_exp_coefs_clim_region))
colnames(neg_exp_coefs_clim_region) <- gsub("\\s", "",colnames(neg_exp_coefs_clim_region))
neg_exp_coefs_clim_region

## power law
pow_coefs_clim_region <- beta.calc.by.ages("pow", "space", "ret", clim=TRUE)
pow_coefs_clim_region$predictor <- gsub("\\d", "", rownames(pow_coefs_clim_region))
colnames(pow_coefs_clim_region) <- gsub("\\s", "",colnames(pow_coefs_clim_region))

# fit measures
neg_exp_aic_clim <- beta.calc.by.ages("neg", "space", "AIC", clim=TRUE)
neg_exp_r2_clim <- beta.calc.by.ages("neg", "space", "pseudo_r2", clim=TRUE)
neg_exp_fit_clim <- data.frame("info"=c("Oligocene", "Miocene", "Pliocene", "Pleistocene", "LastGlacial", "Holocene", "Present"), "pseudo_r2"=neg_exp_r2_clim, "AIC"=neg_exp_aic_clim)
neg_exp_fit_clim

cene_decay_plot <- function(cene, col_, clim=FALSE){
    if (clim == TRUE){
        dist_seq <- seq(from=0, to=30, by=0.2)
        plot_decay <- data.frame("dist"=dist_seq)
        newx <- data.frame("clim_dist" = seq(0, 30, 0.2)) # 下で描画に使う
        icpt <- neg_exp_coefs_clim_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_clim_region$predictor=="(Intercept)"]
        slop <- neg_exp_coefs_clim_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_clim_region$predictor=="clim_dist"]
        plot_decay$pred <- 1-exp((slop * plot_decay$dist + icpt))

        plotdata <- subset(clim_dist_plot_clim_region, clim_dist_plot_clim_region$info==cene) 
        plotdata$dist <- plotdata$clim_dist
        
        plot(beta.sim ~ dist, plotdata, main=paste(cene, ""), xlab="delta temperature [degree]", ylab="turnover", col=col_, pch=20,
        cex.lab  = 2,       #  軸の説明の字の大きさを設定する
        cex.axis = 1.8,      #  軸の数字等（ラベル）の大きさを設定する
        cex.main = 1.8,
        ylim=c(0.0, 0.65), las=1, xlim=c(0, 15))
        
    }else{
        dist_seq <- seq(from=0, to=3000, by=10)
        plot_decay <- data.frame("dist"=dist_seq)
        newx <- data.frame("sp_dist" = seq(0, 3000, 20))
        icpt <- neg_exp_coefs_distonly_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_distonly_region$predictor=="(Intercept)"]
        slop <- neg_exp_coefs_distonly_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_distonly_region$predictor=="sp_dist"]
        plot_decay$pred <- 1-exp((slop * plot_decay$dist + icpt))
        
        plotdata <- subset(plot_data_region_dist, plot_data_region_dist$info==cene) 
        

        plot(beta.sim ~ dist, plotdata, main=paste(cene, ""), xlab="spatial distance [km]", ylab="turnover", col=col_, pch=20,
        cex.lab  = 2,       #  軸の説明の字の大きさを設定する
        cex.axis = 1.8,      #  軸の数字等（ラベル）の大きさを設定する
        cex.main = 1.8,
        ylim=c(0.0, 0.65), las=1, xlim=c(0, 2000))
        
    }



    

    

    if (clim == TRUE){
            m <- beta.regr(cene, plot_data_region_dist, beta, "neg", clim=TRUE)$result   
    }else{
            m <- beta.regr(cene, plot_data_region_dist, beta, "neg")$result
    }

    preds <- predict(m, newdata = newx, se.fit = TRUE, type = "link")
    critval <- 1.96 ## approx 95% CI
    upr <- 1 - exp(preds$fit + (critval * preds$se.fit))
    lwr <- 1 - exp(preds$fit - (critval * preds$se.fit))
    fit <- 1 - exp(preds$fit)
    if (clim == TRUE){
        lines(newx$clim_dist, upr, col = 'darkgreen')
        lines(newx$clim_dist, fit, col = 'orange')
        lines(newx$clim_dist, lwr, col = 'darkgreen')
    }else{
        lines(newx$sp_dist, upr, col = 'darkgreen')
        lines(newx$sp_dist, fit, col = 'orange')
        lines(newx$sp_dist, lwr, col = 'darkgreen')
    }
            
}


options(repr.plot.width=16, repr.plot.height=12)

par(mfrow=c(3,3)) 
for (i in 1:length(geo.ages)){
cene_decay_plot(geo.ages[i], "grey", FALSE)
}


clim_dist_add <- function(cene, data_aggr){
    da <- subset(data_aggr, data_aggr$info == cene)
    if (cene == "LastGlacial"){
        da$clim_dist <- da[,"clim_dist_LGP"]
    }
    else if(cene == "Present"){
        da$clim_dist <- da[,"clim_dist_Present"]
    }
    else{
        da$clim_dist <- da[,paste("clim_dist_", gsub("cene", "", cene), sep="")]        
    }
        return(da)
}


clim_dist_plot_clim_region <- rbind(
clim_dist_add("Oligocene", plot_data_region_dist),
clim_dist_add("Miocene", plot_data_region_dist),
clim_dist_add("Pliocene", plot_data_region_dist),
clim_dist_add("Pleistocene", plot_data_region_dist),
clim_dist_add("LastGlacial", plot_data_region_dist),
clim_dist_add("Holocene", plot_data_region_dist),
clim_dist_add("Present", plot_data_region_dist)
)


options(repr.plot.width=16, repr.plot.height=12)

par(mfrow=c(3,3)) 
for (i in 1:length(geo.ages)){
cene_decay_plot(geo.ages[i], "grey", TRUE)
}







