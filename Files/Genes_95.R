setwd("")
x <- read.table("depth_genes_95_edited.txt", sep="\t", row.names=1, header = TRUE, stringsAsFactors = FALSE)
y <- x[c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27)]

titles <- NULL
titles <- data.frame("Cen01", "Cen03", "Cen04", "Cen05", "Cen06", "Cen07", "Cen10", "Cen12", 
                    "Cen13", "Cen14", "Cen15", "Cen16", "Cen17")
hist <- NULL
for(i in 1:13) {
  hist[[i]] <- hist(y[,i], main = paste("Gene Coverage in", titles[,i]), 
                    xlab = "Gene Coverage")
}

freq <- x[c(21, 25)]
freq
freq[,1]
median_freq <- NULL
for(i in 1:2) {
  median_freq[i] <- median(freq[,i])
}

gene_freq <- freq/median_freq

cen14_freq <- hist(gene_freq[,1], main = "Gene Frequency in Cen14", xlab = "Gene Frequency",
                   col = "indianred1")
cen16_freq <- hist(gene_freq[,2], main = "Gene Frequency in Cen16", xlab = "Gene Frequency",
                   col = "indianred1")

abs_hist <- hist((abs(gene_freq[,1] - gene_freq[,2])), 
                 main="Absolute Value of Gene Frequency Differences \n in Cen14 and Cen16", 
                 xlab = "Magnitude of Gene Frequency Differences",
                 col = "darkmagenta")

non_abs_hist <- hist((gene_freq[,1] - gene_freq[,2]), 
                     main="Gene Frequency Differences \n in Cen14 and Cen16", 
                     xlab = "Gene Frequency Differences",
                     col = "darkmagenta",
                     breaks = 30)
diff_genes <- gene_freq[(abs(gene_freq[,1] - gene_freq[,2])) > 1.0,]
diff_genes
row.names(diff_genes) <- row.names(x[(abs(gene_freq[,1] - gene_freq[,2])) > 1.0,])
diff_genes

