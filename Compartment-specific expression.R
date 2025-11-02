library(ggplot2)
library(ggpubr)
library(ggsci)
library(rstatix)

resolution <- commandArgs(trailingOnly=T)
f <- sprintf("%s.switch.csv", resolution)
out_f <- sprintf("%s.switch.box.pdf", resolution)
df <- read.table(f, sep="\t", header=T)
df <- df[df$compartment_type == "A-to-B" | df$compartment_type == "B-to-A",]
df <- df[df$value < quantile(df$value, 0.90),]
g <- ggboxplot(df, "compartment_type", "value", color = "sample", bxp.errorbar=T, outlier.shape = NA)+
  theme_minimal() + 
  labs(x = "", y = "Gene Expression(TPM)") +  
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
    )  +
  geom_pwc(
    aes(group = sample), 
    method = "wilcox_test", label = "p.adj",
    p.adjust.method = "BH"
  ) 

ggsave(out_f, plot=g, width = 8, height = 6)
