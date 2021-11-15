
# Before beginning, make sure the following packages are installed.
# Required packages
library(zCompositions)
library(coda.base)
library(compositions)
library(ape)
library(vegan)
library(corrplot)
library(fossil)
# You can install packages using BiocManager, example:
#BiocManager::install("vegan")

#Importing Data
#Define community table and metadata map file
# Tell R where to find your metadata file
map_name<-"Moisture_soil.csv" # Change as needed
# Tell R where to find your community data file
table_name<-"Cheatgrass_Bacteria_ESV.csv" # Change input as needed, this example will reference a bacteria dataset

#Import community table
otutablewt <- read.csv(file=table_name, comment.char="", header=T, row.names=1, stringsAsFactors=T, check.names=FALSE)
# check number of taxa and samples in the table
dim(otutablewt)
# Seperate out sample counts from any other metadata in the table
unorder.otutable<-otutablewt[,1:15] # These number reflect the number of samples in the table that I want to seperate from other info.
# Sort the columns by alphanumeric sorting
otutable<-unorder.otutable[,order(colnames(unorder.otutable))]

#Import metadata map file
unorder.map<-read.csv(file=map_name, header=T, comment.char="", , row.names=1, stringsAsFactors=T, check.names=FALSE)
# Sort the rows by alphanumeric sorting
map<-unorder.map[order(rownames(unorder.map)),]
# subset map for only samples found in the community table
sub.map<-subset(map, rownames(map) %in% colnames(otutable))

# Confirm symmetry of metadata and community tables
rownames(sub.map) == colnames(otutable)



#### Betadiv calculations ####
### Bray-curtis distances ###
bray.dist<-vegdist(t(otutable), method="bray") # Accounts for abundances and the presence/absence of each species
### Jaccard distances ###
jac.dist<-vegdist(t(otutable), method="jaccard") # Accounts only for the presence/absence of each species

### Manual Aitchison Method ###
sub.otutable<-t(otutable)
# Trim out ASVs that show up 10x or less
count <- 10 #this is the chosen cutoff
sub.table <- data.frame(sub.otutable[which(apply(sub.otutable, 1, function(x){mean(x)}) > count),],
                                          check.names=F)
# Impute zeros
sub.table_czm <- cmultRepl(sub.table,  label=0, method="CZM")
#sub.table_gbm <- cmultRepl(sub.table,  label=0, method="GBM")                                   
# Log-ratio transform the data
# CLR transformation
P.sub.table_czm_clr <- clr(sub.table_czm)
# Calculate Distances
P.sub.table_czm_clr.dist <- vegdist(P.sub.table_czm_clr, method = "euclidean", diag = FALSE, upper = TRUE)
aitch.dist<-P.sub.table_czm_clr.dist



#### Plot dissmilartiy matrix as NMDS ####
## Define variables for plot code
plot.dm<-metaMDS(t(otutable), distance = "bray")
plot.map<-sub.map[rownames(sub.map) %in% rownames(as.matrix(plot.dm)),] # Change the source of the metadata here
# variance explained by axis (ends up being same as PCoA)
MDS <- cmdscale(vegdist(t(otutable), method = "bray"), k = 2, eig = T, add = T )
round(MDS$eig*100/sum(MDS$eig),1)

pdf("Cheatgrass_bacteria_nmds_bray.pdf") # Name of file output, CHANGE ME as needed
# Plotting 
# FUNCTION:
    Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
      hpts <- chull(x = xcoord, y = ycoord)
      hpts <- c(hpts, hpts[1])
      lines(xcoord[hpts], ycoord[hpts], col = lcolor)
    } 
    # END OF FUNCTION
# Now run the Principle Coordinate Analysis ("nmds") in package ape:
    nmds.pts<-plot.dm$points[,1:2]
    x.minimum<-min(nmds.pts[,1])
    x.maximum<-max(nmds.pts[,1])
    y.minimum<-min(nmds.pts[,2])
    y.maximum<-max(nmds.pts[,2])

    # And to be fancy, b/c Jack put the time into figuring this out,
    # calculate the % variation explained by each axis
    # getting percent explained by each PC vector
    vars_percent <- (nmds.pts$values$Eigenvalues / sum(nmds.pts$values$Eigenvalues)) * 100
# Step 1: make very simple plot, no colors
par(mar = c(5,5,5,1)) 
   
    plot(nmds.pts[,2] ~ nmds.pts[,1],
         xlab=paste("NMDS1: ", round(vars_percent[1], 2), "%", sep=""),
         ylab=paste("NMDS2: ", round(vars_percent[2], 2), "%", sep=""),
         cex=0.1, cex.lab=1.5, cex.axis=1.5, col="white"
         ,main="A. Bacteria", cex.main=2 # CHANGE this to be the title of the plot
         ,xlim=c(x.minimum-0.05, x.maximum+0.2), ylim=c(y.minimum-0.05, y.maximum+0.2)
)
    # Step 2: now subset vectors for the positions in the alpha file that correspond to cohcryo
    map.cheat <- which(plot.map$Treatment == "cheat") ## CHANGE these depending on the number of categories in your variable of interest
    map.native <- which(plot.map$Treatment == "native") ## CHANGE these depending on the number of categories in your variable of interest
   
    # Step 3: and plot the points for cohcryo samples in red over the top using points()
    # Add a set of points and Plot_ConvexHull functions for each category in your variable
    # If may have to remove or add category specific plots as needed
   
    points(nmds.pts[map.cheat,2] ~ nmds.pts[map.cheat,1],
           xlab=paste("nmds1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("nmds2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=19, col="#C00000", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = nmds.pts[map.cheat,1], ycoord = nmds.pts[map.cheat,2], lcolor = "#C00000")

    points(nmds.pts[map.native,2] ~ nmds.pts[map.native,1],
           xlab=paste("nmds1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("nmds2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=19, col="#002060", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = nmds.pts[map.native,1], ycoord = nmds.pts[map.native,2], lcolor = "#002060")

	legend(x.maximum-0.6,y.maximum+0.2, legend=c("Cheatgrass","Native Grasses"), 
		col=c("#C00000","#002060"),
			pch=c(19,19), cex=1.2,
			box.lwd = 0, box.lty=0) # Add labels, colors, and point shape for each category you want to plot
 dev.off()



# Community composition Stats #
am1<-adonis(bray.dist ~ plot.map$Treatment, data = plot.map)
am1
am2<-adonis(aitch.dist ~ plot.map$Treatment, data = plot.map)
am2


#### Alphadiv Calculations ####
# Here are examples of how to get diversity metrics in R
# Shannon Diversity Index
Shannon<-diversity(t(otutable), index = "shannon")
# Simpson Diversity index
Simpson<-diversity(t(otutable), index = "simpson")
# Inverse Simpson (simpson_e) or Evenness Index
Invert.Simpson<-diversity(t(otutable), index = "invsimpson")
# Richness
Richness<-apply(t(otutable)>0,1,sum)
# Chao1
Chao1 <- apply(t(otutable), 1, chao1)

# combine alpha diversity metrics into a single data table
alternative.alphadiv<-as.data.frame(cbind(Shannon,Simpson,Invert.Simpson,Richness,Chao1))
# Confirm symmetry with metadata map data object
rownames(alternative.alphadiv) == rownames(sub.map)

# Export alpha diversity data
write.table(alternative.alphadiv, "bac_alphadiv.txt", sep='\t')

#### Alphadiv Plots ####
# Define plot variables
plot.alpha<-alternative.alphadiv
plot.map<-sub.map

# Plot Box and Whisker plots for alpha-diversity
pdf("Cheatgrass_bacteria_alpha.pdf") # Name of file output, CHANGE ME as needed
for(i in 1:ncol(plot.alpha)){
	name<-colnames(plot.alpha)[i]
	data<-cbind(plot.alpha[,name], plot.map["Treatment"]) # CHANGE this to change which variable in your metadata you want to plot by
	
	if(name == "Richness") {
		par(mar=c(4.5,4.5,3.1,2.1))
		boxplot(data[,1] ~ factor(data[,2]), main="A. Bacteria", ylab="Observed ESVs (per 1g dry soil)", las=1,  # CHANGE title here for plot titles
			col = c("#C00000","#002060"), cex.lab=1.5, cex.main=2, cex.axis=1.1,
			xaxt = "n", xlab="")
			axis(1, at = 1:2, labels = c("Cheatgrass","Native Grasses"), cex.axis = 1.5)
	} else if(name == "PD_whole_tree") {
			par(mar=c(8.1,4.1,4.1,2.1))
		boxplot(data[,1] ~ factor(data[,2]), main=c("A. Bacteria: ",name), ylab="Phylogenetic Diversity", xlab="Site", las=1, # CHANGE title here for plot titles
		#	col = c("#002060","#C00000"))
	} else {
		par(mar=c(8.1,4.1,4.1,2.1))
		boxplot(data[,1] ~ factor(data[,2]), main=c("A. Bacteria: ",name), ylab=paste(name," Index"), xlab="Site", las=1,  # CHANGE title here for plot titles
			col = c("#002060","#C00000"))
	}
}
dev.off()

#### Alphadiversity Stats ###
# t.test
for(i in 1:ncol(plot.alpha)){
	mi<-t.test(plot.alpha[,i] ~ plot.map$Treatment)
	print(colnames(plot.alpha[i]))
	print(summary(mi))
	print(mi)
	print(TukeyHSD(mi))
	}
# Mann - Whitney Test (wilcox.test)
for(i in 1:ncol(plot.alpha)){
print(i)
mi<-wilcox.test(plot.alpha[,i] ~ Treatment, data=plot.map)
print(colnames(plot.alpha[i]))
print(mi)
}


### Extracting Abundant Taxa ###
order.sub.map<-sub.map[order(rownames(sub.map)),]

# Absolute Sequence Abundance
abundsort.otutablewt<-otutablewt[order(otutablewt$Total, decreasing = TRUE),]
sub.abundsort.otutablewt<-as.matrix(abundsort.otutablewt[1:30,order(colnames(abundsort.otutablewt[,1:15]))])
total.tax<-subset(abundsort.otutablewt, select = colnames(abundsort.otutablewt) == "Taxonomy")

# ANOVAs for top 30 abundant taxa
tab<-sub.abundsort.otutablewt
tax<-total.tax

# Calculate species scores and stats from NMDS data object for abundant taxa #
envfit(plot.dm, t(tab))

for(i in 1:nrow(tab)){
# t.test
mi<-t.test(tab[i,] ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(summary(mi))
print(mi)
}
# Mann-Whitney test
for(i in 1:nrow(tab)){
mi<-wilcox.test(tab[i,] ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}


# Relative Abundance
abundsort.otutablewt<-otutablewt[order(otutablewt$Total, decreasing = TRUE),]
prop.mat.otutable<-(prop.table(as.matrix(abundsort.otutablewt[,1:15]), 2)*100)
sub.abundsort.prop.mat.otutable<-prop.mat.otutable[1:30,order(colnames(prop.mat.otutable))]
relative.tax<-subset(abundsort.otutablewt, select = colnames(abundsort.otutablewt) == "taxonomy")

# ANOVAS for Top 30 abundant taxa
tab<-sub.abundsort.prop.mat.otutable
tax<-relative.tax

# Write out relative abundance table to directory
rownames(prop.mat.otutable) == rownames(tax)
relative.table<-cbind(prop.mat.otutable, tax)
#write.table(relative.table, "Cheatgrass_Nematodes_no_Embryophyta_relative_ESV_table.txt", sep='\t')


for(i in 1:nrow(tab)){
# T.test
for(i in 1:nrow(tab)){
mi<-t.test(tab[i,] ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}
# Mann - Whitney Test (wilcox.test)
for(i in 1:nrow(tab)){
mi<-wilcox.test(as.numeric(tab[i,]) ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}

## Significance of Treatment by Phyla ##
#Importing Data
#Define community table and metadata map file
# Tell R where to find you metadata file
map_name<-"Moisture_soil.csv"
# Tell R where to find your community data file
table_name<-"Cheatgrass_Bacteria_Group_ESV.csv"

#Import community table
otutablewt <- read.csv(file=table_name, comment.char="", header=T, row.names=1, stringsAsFactors=T, check.names=FALSE)
dim(otutablewt)
# Seperate out sample counts from any other metadata in the table
unorder.otutable<-otutablewt[,1:15] # These number reflect the number of samples in the table that I want to seperate from other info.
# Sort the columns by alphanumeric sorting
otutable<-unorder.otutable[,order(colnames(unorder.otutable))]

#Import metadata map file
unorder.map<-read.csv(file=map_name, header=T, comment.char="", , row.names=1, stringsAsFactors=T, check.names=FALSE)
# Sort the rows by alphanumeric sorting
map<-unorder.map[order(rownames(unorder.map)),]
# subset map for only samples found in the community table
sub.map<-subset(map, rownames(map) %in% colnames(otutable))

# Confirm symmetry of metadata and community tables
rownames(sub.map) == colnames(otutable)


# Total Abundance
#### Alphadiv  ####
trophic.table<-abundsort.otutablewt[,c(1:15,16)]
trophic.groups<-abundsort.otutablewt[,16]

for(i in 1:length(unique(sort(trophic.groups)))){
group<-unique(sort(trophic.groups))[i]
rows<-which(trophic.table[,16] == group)
sub.table<-trophic.table[rows,1:15]
#print(sub.table)
# Richness
Richness<-apply(sub.table>0,2,sum)
#print(Richness)
# Chao1
Chao1 <- apply(sub.table, 2, chao1)
#print(Chao1)
# Combine metrics into a single data object and order samples alphanumerically
alternative.alphadiv<-as.data.frame(cbind(Richness,Chao1))
order.alternative.alphadiv<-alternative.alphadiv[order(rownames(alternative.alphadiv)),]

# Absolute Abundance 
abundsort.otutablewt<-otutablewt[order(otutablewt$Total, decreasing = TRUE),]
total.trophic<-aggregate(abundsort.otutablewt[,1:15],by=list(abundsort.otutablewt[,16]), FUN=sum) 
order.total.trophic<-total.trophic[,order(colnames(total.trophic))]
final.total.trophic<-order.total.trophic[,1:15]
total.troph<-as.matrix(subset(total.trophic, select = colnames(total.trophic) == "Group.1"))
rownames(final.total.trophic)<-total.troph

# Label new data objects for easier future reference
tab<-final.total.trophic
tax<-total.troph

# Abundance stats
for(i in 1:nrow(tab)){
mi<-aov(as.numeric(tab[i,]) ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(summary(mi))
}
# T.test
for(i in 1:nrow(tab)){
mi<-t.test(as.numeric(tab[i,]) ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}
# Welch Two Sample T.test
for(i in 1:nrow(tab)){
mi<-t.test(as.numeric(tab[i,]) ~ Treatment, data=order.sub.map, alternative = c("two.sided","less","greater"))
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}
# Mann - Whitney Test (wilcox.test)
for(i in 1:nrow(tab)){
mi<-wilcox.test(as.numeric(tab[i,]) ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}


# Relative Abundance
abundsort.otutablewt<-otutablewt[order(otutablewt$Total, decreasing = TRUE),]
aggre<-aggregate(abundsort.otutablewt[,1:15],by=list(abundsort.otutablewt[,16]), FUN=sum) 
prop.mat.otutable<-(prop.table(as.matrix(aggre[,2:16]), 2)*100)
order.prop.mat.otutable<-as.matrix(prop.mat.otutable[,order(colnames(prop.mat.otutable))])
rownames(order.prop.mat.otutable)<-total.troph
tab<-order.prop.mat.otutable
tax<-total.troph

for(i in 1:nrow(tab)){
# T.test
for(i in 1:nrow(tab)){
mi<-t.test(tab[i,] ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}
# Mann - Whitney Test (wilcox.test)
for(i in 1:nrow(tab)){
mi<-wilcox.test(as.numeric(tab[i,]) ~ Treatment, data=order.sub.map)
name<-rownames(tab)[i]
print(rownames(tab)[i])
row<-rownames(tab)[i]
taxonomy<-factor(tax[rownames(tax) == name,])
print(taxonomy)
print(mi)
}


#### Examining Environmental Data ####
### Create a correlation matrix ###
correlations<-cor(sub.map[1:15,17:33], method = "pearson") # These numbers are based on which variables are numeric in our metadata
# Default uses pearson correlations but there are multiple options

# Plot a pairwise depictions of your variable correlations
corrplot(correlations)

# can save as a pdf
pdf("correlations.pdf")
corrplot(correlations)
dev.off()

# PREMANOVA
adonis(formula = bray.dist ~ Percent.Water, data = sub.map) 
adonis(formula = bray.dist ~ Percent.Water + Treatment, data = sub.map) 
adonis(formula = aitch.dist ~ Percent.Water, data = plot.map)
adonis(formula = aitch.dist ~ Percent.Water + Treatment, data = plot.map)
