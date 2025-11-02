library(ggplot2)  
library(ggpubr)
library(ggsci)
library(rstatix)
  
args <- commandArgs(trailingOnly = TRUE)
file <- args[1] #H3K4me3.multi.txt
mark <- args[2]

data <- read.csv(file, stringsAsFactors = FALSE, sep="\t")
data$Value <- log(data$Value + 1)
       
g <- ggboxplot(data, x="Group", y= "Value", fill="comb", outlier.shape=NA, lwd=1, width=0.5) + 
  geom_pwc(
    aes(group = comb), tip.length = 0.01,
    method = "t_test", label = "p.adj.signif",
    bracket.nudge.y = 0,
    step.increase = 0.05,
    vjust=0.5
  ) +
  ggtitle(mark) + 
  theme_minimal() + 
  labs(x = "", y = "Expression level\n(log[TPM values + 1])") +
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
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 16)
    ) + 
    scale_fill_npg() 
  
ggsave(paste0(mark, "_meth_chip_TPM_boxplot.png"), g, width = 6, height = 8, dpi = 300)
