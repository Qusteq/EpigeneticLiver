library(MuMIn)
library(car)
library(relaimpo)

data <- read.table("mat.csv", header=T)
data <- data[data$RNA !=0,]
data <- scale(data, scale = FALSE)
data <- as.data.frame(data)
model <- lm(RNA ~ H3K4me3 + H3K27ac + H3K27me3, data = data)

options(na.action = na.fail)

# select best model using dredge
dr_model <- dredge(model)
summary(dr_model)
best_model <- get.models(dr_model, 1)  
summary(best_model)

# VIF
vif_model <- vif(model)
print(vif_model)

# decompose R^2 using calc.relimp
reimp <- calc.relimp(best_model, type = c("lmg", "last", "first", "pratt"), rela = TRUE)
print(reimp)
