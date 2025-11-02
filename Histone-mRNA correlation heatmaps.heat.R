library(ggplotify)
library(ComplexHeatmap)
library(ggplot2)
library(rcolors)
data <- read.table("mat.csv", header=T)
print(dim(data))
data <- data[data$RNA !=0,]
data <- scale(data, scale = FALSE)
print(dim(data))
mat1 <- as.matrix(data[,4])
colnames(mat1) <- "RNA"
mat2 <- as.matrix(data[,1:3])
colnames(mat2) <- colnames(data)[1:3]


part1 <- get_color(rcolors$CBR_coldhot, 10)[1:5]
part2 <- get_color(rcolors$CBR_coldhot, 20)[11:20]
colors <- c(part1, part2)
svg("5.4.svg", width=5, height=5)
p1 <- pheatmap(mat1, name="RNA-seq NASH\n(log2)", cluster_rows = FALSE, cluster_cols = FALSE, color=colors, angle_col="0") 
p2 <- pheatmap(mat2, name="ChIP-seq NASH\n(log2)", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames=F, color=colors, angle_col="0") 
p1 + p2
dev.off()

