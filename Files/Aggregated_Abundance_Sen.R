setwd("/Users/janelee/Documents/MSU_REU/Coverage_Sensitive/")
library(ggplot2)
library(reshape2)

x <- read.table("CPR_FileNames_Sensitive.txt", stringsAsFactors = FALSE)
y <- NULL
y <- as.list(y)
for (i in 1:12){
  y[[i]] <- read.table(paste(x[i,1]), sep="\t", stringsAsFactors = FALSE, header=TRUE)
}

names(y) <- x[,1]
names(y)
total <- NULL
total <- as.data.frame(total)
for(i in 1:12){
  values <- ((colSums(y[[i]][,seq(4,27,2)]))/nrow(y[[i]]))
  total <- rbind(total, values)
}
total

column_names <- NULL
column_names <- data.frame("Cen01_MA.bam", "Cen03_MA.bam", "Cen04_MA.bam", "Cen05_MA.bam", "Cen06_MA.bam", "Cen07_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam", "Cen15_MA.bam", "Cen16_MA.bam", "Cen17_MA.bam")
column <- t(column_names)
colnames(total) <- column[,1]
row.names(total) <- x[,1]
row.names(total)

total2 <- total[c(1,5,6,7,8,9,10,11,12,2,3,4),]
total2
setwd("/Users/janelee/Documents/MSU_REU/")
metadata <- read.table("Centralia_Collapsed_Map_forR.txt", sep="\t", header = T)
metadata_edit <- metadata[-c(2, 8, 9, 11, 13, 18), ]
metadata_edit

soil_temp <- metadata_edit$SoilTemperature_to10cm
soil_temp

z <- NULL
z <- as.list(z)
for(i in 1:12) {
  z[[i]] <- cor.test(soil_temp, as.numeric(total2[i,]))
}
z

for(i in 1:12) {
  plot(soil_temp, total[i,], xlab = "Temperature °C", ylab = "Abundance")
}

total_transpose <- t(total2)
total_transpose

total_reshape <- melt(total_transpose, id=c("Coverage.METABAT_Sensitive.1", 
  "Coverage.METABAT_Sensitive.2", "Coverage.METABAT_Sensitive.3",
  "Coverage.METABAT_Sensitive.4","Coverage.METABAT_Sensitive.5", "Coverage.METABAT_Sensitive.6",
  "Coverage.METABAT_Sensitive.7", "Coverage.METABAT_Sensitive.8", "Coverage.METABAT_Sensitive.9",
  "Coverage.METABAT_Sensitive.10", "Coverage.METABAT_Sensitive.11", "Coverage.METABAT_Sensitive.12"))
total_reshape

total2
top_4 <- as.data.frame(total_transpose[,1:4])
class(top_4)
class(total_reshape)

p <- ggplot(total_reshape, aes(x = Var1, y= value, group = as.factor(Var2))) 
p + geom_bar(stat = "identity", aes (fill = Var2), position = "dodge") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)
  ) + 
  ggtitle("Coverage in each Centralia Location") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Bins") +
  xlab("Location") + 
  ylab("Coverage") +
  coord_cartesian(ylim=c(0, 20)) +
  scale_fill_discrete(
    breaks=c("Coverage.METABAT_VerySensitive.1", 
             "Coverage.METABAT_VerySensitive.2", "Coverage.METABAT_VerySensitive.3",
             "Coverage.METABAT_VerySensitive.4","Coverage.METABAT_VerySensitive.5", "Coverage.METABAT_VerySensitive.6",
             "Coverage.METABAT_VerySensitive.7", "Coverage.METABAT_VerySensitive.8", "Coverage.METABAT_VerySensitive.9",
             "Coverage.METABAT_VerySensitive.10", "Coverage.METABAT_VerySensitive.11", "Coverage.METABAT_VerySensitive.12"),
    labels=c("Sensitive 1", "Sensitive 2", "Sensitive 3", "Sensitive 4",
             "Sensitive 5", "Sensitive 6", "Sensitive 7", "Sensitive 8",
             "Sensitive 9", "Sensitive 10", "Sensitive 11", "Sensitive 12")) +
  scale_x_discrete(
    breaks=c("Cen01_MA.bam", "Cen03_MA.bam", "Cen04_MA.bam", "Cen05_MA.bam", "Cen06_MA.bam",
             "Cen07_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam",
             "Cen15_MA.bam", "Cen16_MA.bam", "Cen17_MA.bam"),
    labels=c("Cen01", "Cen03", "Cen04", "Cen05", "Cen06",
             "Cen07", "Cen10", "Cen12", "Cen14",
             "Cen15", "Cen16", "Cen17"))



q <- ggplot(total_reshape, aes(x = Var1, y= value, group = as.factor(Var2))) 
q + geom_bar(stat = "identity", aes (fill = Var2), position = "dodge") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)
  ) + 
  ggtitle("Coverage in each Centralia Location") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Bins") +
  xlab("Location") + 
  ylab("Coverage") +
  coord_cartesian(ylim=c(0, 0.03)) +
  scale_fill_discrete(
    breaks=c("Coverage.METABAT_VerySensitive.1", 
             "Coverage.METABAT_VerySensitive.2", "Coverage.METABAT_VerySensitive.3",
             "Coverage.METABAT_VerySensitive.4","Coverage.METABAT_VerySensitive.5", "Coverage.METABAT_VerySensitive.6",
             "Coverage.METABAT_VerySensitive.7", "Coverage.METABAT_VerySensitive.8", "Coverage.METABAT_VerySensitive.9",
             "Coverage.METABAT_VerySensitive.10", "Coverage.METABAT_VerySensitive.11", "Coverage.METABAT_VerySensitive.12"),
    labels=c("Sensitive 1", "Sensitive 2", "Sensitive 3", "Sensitive 4",
             "Sensitive 5", "Sensitive 6", "Sensitive 7", "Sensitive 8",
             "Sensitive 9", "Sensitive 10", "Sensitive 11", "Sensitive 12"))  +
  scale_x_discrete(
    breaks=c("Cen01_MA.bam", "Cen03_MA.bam", "Cen04_MA.bam", "Cen05_MA.bam", "Cen06_MA.bam",
             "Cen07_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam",
             "Cen15_MA.bam", "Cen16_MA.bam", "Cen17_MA.bam"),
    labels=c("Cen01", "Cen03", "Cen04", "Cen05", "Cen06",
             "Cen07", "Cen10", "Cen12", "Cen14",
             "Cen15", "Cen16", "Cen17"))

soil_df <- as.data.frame(soil_temp)
soil_df
df <- data.frame(x=1:12, y=soil_df$soil_temp)
df
soil_plot <- ggplot(df, aes(x,y))
soil_plot + geom_line(stat = "identity", aes (x,y)) + 
  scale_x_discrete(
    labels=c("Cen01", "Cen03", "Cen04", "Cen05", "Cen06",
             "Cen07", "Cen10", "Cen11", "Cen12", "Cen14",
             "Cen15", "Cen16", "Cen17"))

soil_plot <- ggplot(df, aes(x,y)) 
soil_plot + geom_line(stat = "identity", aes (x,y)) +
  scale_x_continuous(breaks=1:12, labels=c("Cen01", "Cen03", "Cen04", "Cen05", "Cen06",
                                           "Cen07", "Cen10", "Cen12", "Cen14","Cen15", "Cen16", "Cen17")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  xlab("Location") + 
  ylab("Temperature °C") +
  ggtitle("Soil Temperatures in each Centralia Location") +
  theme(plot.title = element_text(hjust = 0.5)) 


total_reshape
a <- ggplot(subset(total_reshape, Var2 %in% c("Coverage.METABAT_VerySensitive.1", "Coverage.METABAT_VerySensitive.2",
      "Coverage.METABAT_VerySensitive.3", "Coverage.METABAT_VerySensitive.4")), 
      aes(x = Var1, y= value, group = as.factor(Var2))) 
a + geom_bar(stat = "identity", aes (fill = Var2), position = "dodge") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)
  ) + 
  ggtitle("Coverage in each Centralia Location: Top 4 Sensitive Bins") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Bins") +
  xlab("Location") + 
  ylab("Coverage") +
  coord_cartesian(ylim=c(0, 20)) +
  scale_fill_discrete(
    breaks=c("Coverage.METABAT_VerySensitive.1", 
             "Coverage.METABAT_VerySensitive.2", "Coverage.METABAT_VerySensitive.3",
             "Coverage.METABAT_VerySensitive.4"),
    labels=c("Sensitive 1", "Sensitive 2", "Sensitive 3", "Sensitive 4")) +
  scale_x_discrete(
    breaks=c("Cen01_MA.bam", "Cen03_MA.bam", "Cen04_MA.bam", "Cen05_MA.bam", "Cen06_MA.bam",
             "Cen07_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam",
             "Cen15_MA.bam", "Cen16_MA.bam", "Cen17_MA.bam"),
    labels=c("Cen01", "Cen03", "Cen04", "Cen05", "Cen06",
             "Cen07", "Cen10", "Cen12", "Cen14",
             "Cen15", "Cen16", "Cen17"))

total_reshape
b <- ggplot(subset(total_reshape, Var1 %in% c("Cen01_MA.bam","Cen03_MA.bam", "Cen04_MA.bam", 
  "Cen05_MA.bam", "Cen07_MA.bam", "Cen17_MA.bam")), aes(x = Var1, y= value, group = as.factor(Var2))) 
b + geom_bar(stat = "identity", aes (fill = Var2), position = "dodge") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)
  ) + 
  ggtitle("Coverage in Cold Sites") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Bins") +
  xlab("Location") + 
  ylab("Coverage") +
  coord_cartesian(ylim=c(0, 0.03)) +
  scale_x_discrete(
    breaks=c("Cen01_MA.bam", "Cen03_MA.bam", "Cen04_MA.bam", "Cen05_MA.bam",
             "Cen07_MA.bam", "Cen17_MA.bam"),
    labels=c("Cen01", "Cen03", "Cen04", "Cen05",
             "Cen07", "Cen17")) +
  scale_fill_discrete(
    breaks=c("Coverage.METABAT_VerySensitive.1", 
             "Coverage.METABAT_VerySensitive.2", "Coverage.METABAT_VerySensitive.3",
             "Coverage.METABAT_VerySensitive.4","Coverage.METABAT_VerySensitive.5", "Coverage.METABAT_VerySensitive.6",
             "Coverage.METABAT_VerySensitive.7", "Coverage.METABAT_VerySensitive.8", "Coverage.METABAT_VerySensitive.9",
             "Coverage.METABAT_VerySensitive.10", "Coverage.METABAT_VerySensitive.11", "Coverage.METABAT_VerySensitive.12"),
    labels=c("Sensitive 1", "Sensitive 2", "Sensitive 3", "Sensitive 4",
             "Sensitive 5", "Sensitive 6", "Sensitive 7", "Sensitive 8",
             "Sensitive 9", "Sensitive 10", "Sensitive 11", "Sensitive 12"))

c <- ggplot(subset(total_reshape, Var1 %in% c("Cen06_MA.bam","Cen10_MA.bam", "Cen12_MA.bam", 
  "Cen14_MA.bam", "Cen15_MA.bam", "Cen16_MA.bam")), aes(x = Var1, y= value, group = as.factor(Var2))) 
c + geom_bar(stat = "identity", aes (fill = Var2), position = "dodge") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)
  ) + 
  ggtitle("Coverage in Hot Sites") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Bins") +
  xlab("Location") + 
  ylab("Coverage") +
  coord_cartesian(ylim=c(0, 0.03)) +
  scale_x_discrete(
    breaks=c("Cen06_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam",
             "Cen15_MA.bam", "Cen16_MA.bam"),
    labels=c("Cen06", "Cen10", "Cen12", "Cen14",
             "Cen15", "Cen16")) +
  scale_fill_discrete(
    breaks=c("Coverage.METABAT_VerySensitive.1", 
             "Coverage.METABAT_VerySensitive.2", "Coverage.METABAT_VerySensitive.3",
             "Coverage.METABAT_VerySensitive.4","Coverage.METABAT_VerySensitive.5", "Coverage.METABAT_VerySensitive.6",
             "Coverage.METABAT_VerySensitive.7", "Coverage.METABAT_VerySensitive.8", "Coverage.METABAT_VerySensitive.9",
             "Coverage.METABAT_VerySensitive.10", "Coverage.METABAT_VerySensitive.11", "Coverage.METABAT_VerySensitive.12"),
    labels=c("Sensitive 1", "Sensitive 2", "Sensitive 3", "Sensitive 4",
             "Sensitive 5", "Sensitive 6", "Sensitive 7", "Sensitive 8",
             "Sensitive 9", "Sensitive 10", "Sensitive 11", "Sensitive 12"))
 

soil_temp
