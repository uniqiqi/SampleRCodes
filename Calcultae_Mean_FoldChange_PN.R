setwd("~/Documents/Qiqi-2022-work/2022-04/#2 BGC")
BGC_filename <- "BGC_class_matrix_PN.csv"
genome_filename <- "ALL_database_0428_PN.csv"

levels = c("Phylum", "Class", "Order", "Family", "Genus")
p = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria")
o = c("Corynebacteriales", "Micrococcales", "Micromonosporales", "Streptomycetales",
      "Bacillales", "Lactobacillales", "Rhizobiales", "Rhodospirillales",
      "Sphingomonadales", "Burkholderiales", "Enterobacterales", "Pseudomonadales")
g = c("Flavobacterium", "Bacillus")
f = c("Nocardiaceae", "Microbacteriaceae","Micrococcaceae", "Flavobacteriaceae",
      "Weeksellaceae", "Bacillaceae","Enterobacteriaceae", "Yersiniaceae")
c = c("Flavobacteriia", "Bacilli","Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")

OTU = list(p,c,o,f,g)

DT <- NULL
COLNM <- NULL
n = 0

for(l in 1:5){
  level = levels[l]
  for(otu in OTU[[l]]){
    n <- n+1
    dt <- NULL
    #find # pa and npa genomes in the designated otu
    genome <- read.csv(genome_filename, header = T)
    genome <- genome[which(genome[,level] == otu),]
    pa <- genome[which(genome$GOLD.Ecosystem.Category == "Plants"),]
    npa <- genome[which(genome$GOLD.Ecosystem.Category != "Plants"),]
    pa <- length(pa[,1])
    npa <- length(npa[,1])
    
    #find # pa and npa BGCs in the designated otu
    genome <- read.csv(genome_filename, header = T)
    BGC <- read.csv(BGC_filename, header = T, row.names = 1)
    genome <- genome[match(rownames(BGC),genome$filename),] # remove genome without any BGC
    genome <- genome[which(genome[,level] == otu),]
    genome <- droplevels(genome)
    BGC <- BGC[match(genome$filename,rownames(BGC)),]
    # result: both genome metadata and BGC matrix contain only genomes from designated otu
    # remove BGC classes that have zero prediction in the current otu
    BGC <- BGC[,which(colSums(BGC)!=0)]
    genome <- data.frame(row.names = genome$filename, trait = genome$GOLD.Ecosystem.Category)
    BGC <- merge(BGC, genome, by = 0)
    PA_BGC <- BGC[which(BGC$trait == "Plants"),]
    PA_BGC[,1] <- NULL
    PA_BGC[,length(PA_BGC)] <- NULL
    NPA_BGC <- BGC[which(BGC$trait != "Plants"),]
    NPA_BGC[,1] <- NULL
    NPA_BGC[,length(NPA_BGC)] <- NULL
    
    dt <- data.frame(colSums(PA_BGC), colSums(NPA_BGC), ((colSums(PA_BGC)/pa)/(colSums(NPA_BGC)/npa)))
    colnm <- c(paste(otu, "_PA", sep = ""), paste(otu, "_NPA", sep = ""), paste(otu, "_PA/NPA", sep = ""))
    colnames(dt) <- colnm
    if(n == 1){
      DT <- dt
      COLNM <- colnm
    }else{
      DT <- merge(DT, dt, by = 0, all = T)
      row.names(DT) <- DT[,1]
      DT[,1] <- NULL
    }
  }
}
#write.csv(DT, "BGC_Mean_Fold_Change_PA.csv", row.names = F, quote = F)
d <- select(DT, ends_with("_PA/NPA"))
write.csv(d, "BGC class enrichment/HPG/BGC_Mean_Fold_Change_PN.csv",quote = F)
  
