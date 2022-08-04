library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
Sys.setlocale(category = "LC_ALL", locale = "Greek")

setwd("~/Documents/Qiqi-2022-work/2022-04/#3 Class-enrichment/HPG")
data <- read.csv("PN_HPG_heatmap.csv", header = T, row.names = 1)
#data[is.na(data)] <- 0
m <- data[,c(9:45)]
i <- data[,c(8,1:7)]
w = length(m[1,])
h = length(m[,1])

### change numbers to levels
#m[m >=5] <- ">5"
#m[m < 5 & m >= 2.5] <- "2.5-5"
#m[m < 2.5 & m >= 2] <- "2-2.5"
#m[m < 2 & m >= 1.5] <- "1.5-2"
#m[m < 1.5 & m >= 1.05] <- "1.05-1.5"
#m[m < 1.05 & m > 1] <- "Non-significant/Low Fold Change/Not Detected"
#m[m == 0] <- "Non-significant/Low Fold Change/Not Detected"
#m[m <= 0.2 & m > 0] <- "0.67-0.95"
#m[m < 0.4 & m >= 0.2] <- "0.5-0.67"
#m[m < 0.5 & m >= 0.4] <- "0.4-0.5"
#m[m < 0.67 & m >= 0.5] <- "0.2-0.4"
#m[m < 0.95 & m > 0.67] <- "<0.2"
#m[m < 1 & m > 0.95] <- "Non-significant/Low Fold Change/Not Detected"

colors = c("#006633","#41874d","#6fa969","#9dcc87","#cbf0a7", "white",
                "#a7daf0","#78c0ea","#4fa4e3","#3286db","#3366cc")
levels = c(">5", "2.5-5", "2-2.5", "1.5-2", "1.05-1.5",
           "Non-significant/Low Fold Change/Not Detected",
          "0.67-0.95","0.5-0.67", "0.4-0.5", "0.2-0.4", "<0.2")
col_fun = colorRamp2(c(5,2.5,2,1.5,1.05,1,0.95,0.67,0.5,0.4,0.2), colors)
lgd = Legend(labels = levels, legend_gp = gpar(fill = colors), title = "Fold Change of Average BGC per Genome")


###Taxa group
groups = c("A", "B", "C", "D")
samples = c(8,5,6,12)
taxa_groups <- NULL
for(n in seq(1, length(groups))){
  group_n = groups[n]
  sample_n = samples[n]
  g = c(rep(group_n, sample_n))
  taxa_groups = c(taxa_groups, g)}
taxa_cols <- c("#EF476F", "#FFD166", "#06D6A0", "#118AB2")

###Taxa annotaion
taxa_info <- i$levels
ha = rowAnnotation(foo = anno_text(
  taxa_info,gp = gpar(col = rep(taxa_cols, samples), 
                      fill = "white", border = FALSE, fontsize = 9)))
ha2 = rowAnnotation('Analyzed BGC' = anno_barplot(cbind(i$P.b, i$N.b), width = unit(4, "cm"), bar_width = 0.8,
                                          gp = gpar(fill = c("#9dcc87","#78c0ea" ) , col = c("#9dcc87","#78c0ea"))))
lgd2 = Legend(labels = c("Plant-associated BIGs", "Non-plant-associated BIGs"), 
             title = " ", legend_gp = gpar(fill = c("#9dcc87","#78c0ea" )))
heatmap = Heatmap(m, name = "Fold change of average BGC/genome", col = col_fun, na_col = "white",
        cluster_rows = F, cluster_columns = F,
        row_names_gp = gpar(fontsize = 0),
        column_names_side = "top", column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
        row_split = taxa_groups, row_title = NULL,
        rect_gp = gpar(col = "black", lwd = 1),
        width = unit(w/2.5, "cm"), height = unit(h/2.5, "cm"),
        left_annotation = ha, right_annotation = ha2, show_heatmap_legend = FALSE)

svg(filename = "Class_enrichment_HPG_PN2.svg", 
    width = unit(w+5, "cm"), 
    height = unit(h+5, "cm"))
draw(heatmap)
draw(lgd, x = unit(74, "cm"), y = unit(48, "cm"))
draw(lgd2, x = unit(72.4, "cm"), y = unit(45, "cm"))
dev.off()







