library(ggplot2)
library(reshape2)

setwd("")

x <- read.table("CPR_FileNames_Specific.txt", stringsAsFactors = FALSE)
y <- NULL
y <- as.list(y)
for (i in 1:4){
  y[[i]] <- read.table(paste(x[i,1]), sep="\t", stringsAsFactors = FALSE, header=TRUE)
}

names(y) <- x[,1]
names(y)

total <- NULL
total <- as.data.frame(total)
for(i in 1:4){
  values <- ((colSums(y[[i]][,seq(4,27,2)]))/nrow(y[[i]]))
  total <- rbind(total, values)
}

column_names <- NULL
column_names <- data.frame("Cen01_MA.bam", "Cen03_MA.bam", "Cen04_MA.bam", "Cen05_MA.bam", "Cen06_MA.bam", "Cen07_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam", "Cen15_MA.bam", "Cen16_MA.bam", "Cen17_MA.bam")
column <- t(column_names)
colnames(total) <- column[,1]

row.names(total) <- x[,1]

metadata <- read.table("Centralia_Collapsed_Map_forR.txt", sep="\t", header = T)
metadata_edit <- metadata[-c(2, 8, 9, 11, 13, 18), ]

soil_temp <- metadata_edit$SoilTemperature_to10cm

z <- NULL
z <- as.list(z)
for(i in 1:4) {
  z[[i]] <- cor.test(soil_temp, as.numeric(total[i,]))
}

for(i in 1:4) {
plot(soil_temp, total[i,], xlab = "Temperature °C", ylab = "Abundance")
}

total_transpose <- t(total)

total_reshape <- melt(total_transpose, id=c("Coverage.METABAT_VerySpecific_Trial.1", 
  "Coverage.METABAT_VerySpecific_Trial.2", "Coverage.METABAT_VerySpecific_Trial.3",
  "Coverage.METABAT_VerySpecific_Trial.4"))

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
    breaks=c("Coverage.METABAT_VerySpecific_Trial.1", "Coverage.METABAT_VerySpecific_Trial.2", 
             "Coverage.METABAT_VerySpecific_Trial.3", "Coverage.METABAT_VerySpecific_Trial.4"),
    labels=c("Specific 1", "Specific 2", "Specific 3", "Specific 4")) +
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
    breaks=c("Coverage.METABAT_VerySpecific_Trial.1", "Coverage.METABAT_VerySpecific_Trial.2", 
             "Coverage.METABAT_VerySpecific_Trial.3", "Coverage.METABAT_VerySpecific_Trial.4"),
    labels=c("Specific 1", "Specific 2", "Specific 3", "Specific 4")) +
  scale_x_discrete(
    breaks=c("Cen01_MA.bam", "Cen03_MA.bam", "Cen04_MA.bam", "Cen05_MA.bam", "Cen06_MA.bam",
             "Cen07_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam",
             "Cen15_MA.bam", "Cen16_MA.bam", "Cen17_MA.bam"),
    labels=c("Cen01", "Cen03", "Cen04", "Cen05", "Cen06",
             "Cen07", "Cen10", "Cen12", "Cen14",
             "Cen15", "Cen16", "Cen17"))

soil_df <- as.data.frame(soil_temp)
df <- data.frame(x=1:12, y=soil_df$soil_temp)
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

cold_sites <- ggplot(subset(total_reshape, Var1 %in% c("Cen01_MA.bam","Cen03_MA.bam", "Cen04_MA.bam", 
                                              "Cen05_MA.bam", "Cen07_MA.bam", "Cen17_MA.bam")), aes(x = Var1, y= value, group = as.factor(Var2))) 
cold_sites + geom_bar(stat = "identity", aes (fill = Var2), position = "dodge") +
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
    breaks=c("Coverage.METABAT_VerySpecific_Trial.1", 
             "Coverage.METABAT_VerySpecific_Trial.2", "Coverage.METABAT_VerySpecific_Trial.3",
             "Coverage.METABAT_VerySpecific_Trial.4"),
    labels=c("Specific 1", "Specific 2", "Specific 3", "Specific 4"))

hot_sites <- ggplot(subset(total_reshape, Var1 %in% c("Cen06_MA.bam","Cen10_MA.bam", "Cen12_MA.bam", 
     "Cen14_MA.bam", "Cen15_MA.bam", "Cen16_MA.bam")), aes(x = Var1, y= value, group = as.factor(Var2))) 
hot_sites + geom_bar(stat = "identity", aes (fill = Var2), position = "dodge") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)
  ) + 
  ggtitle("Coverage in Hot Sites") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Bins") +
  xlab("Location") + 
  ylab("Coverage") +
  coord_cartesian(ylim=c(0, 0.03)) +
  scale_fill_discrete(
    breaks=c("Coverage.METABAT_VerySpecific_Trial.1", 
             "Coverage.METABAT_VerySpecific_Trial.2", "Coverage.METABAT_VerySpecific_Trial.3",
             "Coverage.METABAT_VerySpecific_Trial.4"),
    labels=c("Specific 1", "Specific 2", "Specific 3", "Specific 4")) +
  scale_x_discrete(
    breaks=c("Cen06_MA.bam", "Cen10_MA.bam", "Cen12_MA.bam", "Cen14_MA.bam",
           "Cen15_MA.bam", "Cen16_MA.bam"),
    labels=c("Cen06", "Cen10", "Cen12", "Cen14", "Cen15", "Cen16")) 
