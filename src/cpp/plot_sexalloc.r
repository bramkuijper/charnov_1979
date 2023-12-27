#!/usr/bin/env Rscript

library("tidyverse")
library("patchwork")

the_data <- read_csv(file="data.csv")

p1 <- ggplot(data=the_data
        ,mapping=aes(x=time, y=meanr)) +
    geom_line() +
    ylim(0,1) +
    theme_classic()

p2 <- ggplot(data=the_data
        ,mapping=aes(x=time, y=varr)) +
    geom_line() +
    theme_classic()

the_data_l <- the_data %>% select(c(time,meanw_gain_curvem,meanw_gain_curvef)) %>%
pivot_longer(cols=c(meanw_gain_curvem,meanw_gain_curvef),names_to="fitness_type"
        ,values_to="fitness")

p3 <- ggplot(data=the_data_l
        ,mapping=aes(x=time, y=fitness)) +
    geom_line(mapping=aes(colour=fitness_type)) +
    theme_classic()

the_data_l <- the_data %>% select(c(time,meanwf,meanwm)) %>%
    pivot_longer(cols=c(meanwf,meanwm),names_to="fitness_type"
            ,values_to="fitness")

p4 <- ggplot(data=the_data_l
        ,mapping=aes(x=time, y=fitness)) +
    geom_line(mapping=aes(colour=fitness_type)) +
    theme_classic()

p1 / p2  / p3 / p4

ggsave(filename="rplot.pdf")
