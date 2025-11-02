library(jsonlite)
library(UpSetR)
library(ggplotify)
library(ggplot2)

json_file <- 'upset_data.json'
data <- fromJSON(json_file)
p <- as.ggplot(upset(fromList(data), order.by = "freq"))
ggsave("upset.svg", p) 
