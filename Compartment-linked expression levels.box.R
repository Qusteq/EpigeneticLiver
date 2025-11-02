library(ggplot2)  
library(ggpubr)
library(ggsci)
library(rstatix)
  
args <- commandArgs(trailingOnly = TRUE)
file <- args[1] #H3K4me3.multi.txt
out_fig <- args[2]

data <- read.csv(file, stringsAsFactors = FALSE, sep="\t")
data$value <- log(data$value + 1)
g <- ggboxplot(data, x="compartment_type", y= "value", fill="compartment_type", outlier.shape=NA, lwd=1, width=0.5, order = c("Total", "A", "B")) + 
  geom_pwc(
    aes(group = compartment_type), tip.length = 0,
    method = "wilcox_test", label="p.adj.format", #label = "p.adj.signif",
    p.adjust.method = "bonferroni", p.adjust.by = "panel"
  ) +
  theme_minimal() + 
  labs(x = "", y = "log(Average TPM)") + 
  border(size = 2) +
  theme(
    plot.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
	legend.position="none"
    ) + 
    scale_fill_npg() 
  
ggsave(out_fig, g, width = 6, height = 8, dpi = 300)
