#---Initialization-----
library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)

isosegmenter <- read.csv("isosegmenter_scoring_error5.csv", header = TRUE)
isoplotter <- read.csv("isoplotter_scoring_error5.csv", header = TRUE)
eq_pred_lengths <- read.csv("EqDom_length_correl.csv", header = TRUE)
var_pred_lengths <- read.csv("VarDom_length_correl.csv", header = TRUE)

print("sensitivity test")
print(wilcox.test(isoplotter$Sens, isosegmenter$Sens))
print("PPV/precsion test")
print(wilcox.test(isoplotter$PPV, isosegmenter$PPV))
print("Jaccard test")
print(wilcox.test(isoplotter$Jaccard, isosegmenter$Jaccard))


isoseg.eq <-isosegmenter %>% 
            filter(TYPE == "E") %>% 
            select(DOMAINS, Sens, PPV, Jaccard) %>% 
            group_by(DOMAINS) %>% 
            summarize(mean_s = mean(Sens),
                      mean_p = mean(PPV),
                      mean_j = mean(Jaccard))

isoseg.var <-isosegmenter %>% 
  filter(TYPE == "V") %>% 
  select(DOMAINS, Sens, PPV, Jaccard) %>% 
  group_by(DOMAINS) %>% 
  summarize(mean_s = mean(Sens),
            mean_p = mean(PPV),
            mean_j = mean(Jaccard))

isoplotter.eq <-isoplotter %>% 
  filter(TYPE == "E") %>% 
  select(DOMAINS, Sens, PPV, Jaccard) %>% 
  group_by(DOMAINS) %>% 
  summarize(mean_s = mean(Sens),
            mean_p = mean(PPV),
            mean_j = mean(Jaccard))

isoplotter.var <-isoplotter %>% 
  filter(TYPE == "V") %>% 
  select(DOMAINS, Sens, PPV, Jaccard) %>% 
  group_by(DOMAINS) %>% 
  summarize(mean_s = mean(Sens),
            mean_p = mean(PPV),
            mean_j = mean(Jaccard))

common_theme <- function() {
  theme_bw() + 
  theme(plot.title=element_text(size=14, lineheight=0.8, hjust=0.5),
        legend.justification = c(1,0))
}

#----Performance Plots-----

plot.isoseg.eq <- ggplot(data = isoseg.eq, aes(x = DOMAINS)) +
    geom_line(aes(y = mean_s, color = "Sensitivity"), size = 2) +
    geom_line(aes(y = mean_p, color = "Precision"), size = 2) +
    geom_line(aes(y = mean_j, color = "Jaccard"), size = 2) +
    labs(x = "Number of Domains", y = "Performance Metrics") +
    ggtitle("(A) isoSegmenter: Equal-Length Domains") + 
    scale_color_manual("Metric", 
                     breaks = c("Sensitivity", "Precision", "Jaccard"),
                     values = c("blue", "red", "green")) + 
    common_theme()+
    theme(legend.position = c(0.95,0.3))

plot.isoseg.var <- ggplot(data = isoseg.var, aes(x = DOMAINS)) +
  geom_line(aes(y = mean_s, color = "Sensitivity"), size = 2) +
  geom_line(aes(y = mean_p, color = "Precision"), size = 2) +
  geom_line(aes(y = mean_j, color = "Jaccard"), size = 2) +
  labs(x = "Number of Domains", y = "Performance Metrics") +
  ggtitle("(B) isoSegmenter: Variable-Length Domains") + 
  scale_color_manual("Metric",
                     breaks = c("Sensitivity", "Precision", "Jaccard"),
                     values = c("blue", "red", "green")) + 
  common_theme()+
  theme(legend.position = c(0.95,0.3))


plot.isoplotter.eq <- ggplot(data = isoplotter.eq, aes(x = DOMAINS)) +
  geom_line(aes(y = mean_s, color = "Sensitivity"), size = 2) +
  geom_line(aes(y = mean_p, color = "Precision"), size = 2) +
  geom_line(aes(y = mean_j, color = "Jaccard"), size = 2) +
  labs(x = "Number of Domains", y = "Performance Metrics") +
  ggtitle("(A) isoPlotter: Equal-Length Domains") + 
  scale_color_manual("Metric", 
                    breaks = c("Sensitivity", "Precision", "Jaccard"),
                    values = c("blue", "red", "green")) + 
  common_theme()+
  theme(legend.position = c(0.95,0))

plot.isoplotter.var <- ggplot(data = isoplotter.var, aes(x = DOMAINS)) +
  geom_line(aes(y = mean_s, color = "Sensitivity"), size = 2) +
  geom_line(aes(y = mean_p, color = "Precision"), size = 2) +
  geom_line(aes(y = mean_j, color = "Jaccard"), size = 2) +
  labs(x = "Number of Domains", y = "Performance Metrics", size = 2) +
  ggtitle("(B) isoPlotter: Variable-Length Domains") +
  scale_color_manual("Metric", 
                     breaks = c("Sensitivity", "Precision", "Jaccard"),
                     values = c("blue", "red", "green")) + 
  common_theme()+
  theme(legend.position = c(0.95,0.3))



#----Precision and Sensitivity Plot-------

one <- data.frame(algorithm = "isoPlotter", mean_s = isoplotter.eq$mean_s, mean_p = isoplotter.eq$mean_p, DOMAINS = isoplotter.eq$DOMAINS)
two <- data.frame(algorithm = "isoSegmenter", mean_s = isoseg.eq$mean_s, mean_p = isoseg.eq$mean_p, DOMAINS = isoseg.eq$DOMAINS)
scatter.data <- rbind(one, two)

scatter <- ggplot(data = scatter.data, aes(x = mean_p, y = mean_s, color = algorithm, label = DOMAINS)) +
  geom_point() +
  geom_text_repel(data=scatter.data, max.overlaps = 20) +
  geom_segment(aes(xend = c(tail(mean_p, n=-1), NA), yend = c(tail(mean_s, n=-1), NA))) +
  labs(x = "Precision", y = "Sensitivity") +
  theme_bw() +
  scale_color_manual(values = c("blue", "red"))


#---------Ground truth to prediction stuff------
require(MASS)
require(scales)
data(Animals)

plot.lengthvar.eq <- ggplot(data = eq_pred_lengths, aes(x = TRUE_LENGTH)) +
  geom_point(aes(y = ISOSEGMENTER_AVG, color = "isoSegmenter"), size = 3) +
  geom_point(aes(y = ISOPLOTTER_AVG, color = "isoPlotter"), size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "True Domain Length", y = "Average Predicted Domain Length", size = 2) +
  ggtitle("Average Predicted Domain Length versus Ground Truth (equal-length)") +
  scale_color_manual("Legend", 
                     breaks = c("isoSegmenter", "isoPlotter"),
                     values = c("blue", "red")) + 
  common_theme() +
  annotation_logticks()

plot.numvar.var <- ggplot(data = var_pred_lengths, aes(x = TRUE_DOMAINS)) +
  geom_point(aes(y = ISOSEGMENTER_NUM, color = "isoSegmenter"), size = 3) +
  geom_point(aes(y = ISOPLOTTER_NUM, color = "isoPlotter"), size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Number of Domains", y = "Avg. Predicted Number of Domains", size = 2) +
  ggtitle("Average Predicted Number of Domains versus Ground Truth (variable-length)") +
  scale_color_manual("Legend", 
                     breaks = c("isoSegmenter", "isoPlotter"),
                     values = c("blue", "red")) + 
  common_theme()

#----Distribution of predictions----

a <- data.frame(algorithm = "isoSegmenter", value = var_pred_lengths$ISOSEGMENTER_AVG)
b <- data.frame(algorithm = "isoPlotter", value = var_pred_lengths$ISOPLOTTER_AVG)
dist.data.var <- rbind(a, b)

plot.lengthdist.var <- ggplot(data = dist.data.var, aes(x = value, fill=algorithm)) +
  geom_density(color = "black", alpha=0.4) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Avg. Predicted Domain Lengths", size = 2) +
  ggtitle("(B) Variable-length Predictions") +
  theme_bw()+
  annotation_logticks(sides="b")+
  theme(legend.justification = c(0,1), legend.position = c(0.05,0.95))

c <- data.frame(algorithm = "isoSegmenter", value = eq_pred_lengths$ISOSEGMENTER_AVG)
d <- data.frame(algorithm = "isoPlotter", value = eq_pred_lengths$ISOPLOTTER_AVG)
dist.data.eq <- rbind(c, d)

plot.lengthdist.eq <- ggplot(data = dist.data.eq, aes(x = value, fill=algorithm)) +
  geom_density(color="black", alpha=0.4) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Avg. Predicted Domain Lengths (bp)", size = 2) +
  ggtitle("(A) Equal-length Predictions") +
  theme_bw()+
  annotation_logticks(sides="b")+
  theme(legend.justification = c(0,1), legend.position = c(0.05,0.95))


#-----Print plots----
grid.arrange(plot.isoseg.eq, plot.isoseg.var, ncol=2)
grid.arrange(plot.isoplotter.eq, plot.isoplotter.var, ncol=2)
print(plot.lengthvar.eq)
grid.arrange(plot.lengthdist.eq, plot.lengthdist.var, ncol=2)