library(networkD3)
library(webshot)

nodes <- read.table("nodes.txt", sep="\t", header=T)
links <- read.table("edges.txt", sep="\t", header=T)

p <- sankeyNetwork(Links = links, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "value", NodeID = "Node",
 fontSize= 12, nodeWidth = 30, sinksRight=FALSE)

saveNetwork(p, "sankey.html")
webshot("sankey.html" , "sankey.svg")
