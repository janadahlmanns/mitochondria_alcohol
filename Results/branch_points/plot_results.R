library(ggplot2)
library(tidyverse)

data <- read.csv("./Analysis/branch_points/counts.csv", sep=";")

data %>% ggplot(aes(x = condition, y = ratio.of.mito.positive.branch.points)) + geom_boxplot(aes(group = condition)) + geom_point()

