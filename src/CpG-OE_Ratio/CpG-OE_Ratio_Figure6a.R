# The script wes used to plot the violine figure.

# Example data:
# tissue  degree  score
# liver   m6a     0.446706
# liver   m6a     0.199786
# liver   m6a     0.406499
# liver   m6a     0.824236
# liver   m6a     0.908658
# liver   m6a     0.649302

setwd("")
input_file <- ""
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
#
data <- read.table(input_file, sep = "\t", header = FALSE)
names(data) <- c("tissue", "degree", "score")
head(data)
#
data$tissue <- capFirst(data$tissue)
data$degree <- factor(data$degree, levels = c("free", "m6a"), ordered = TRUE)
data$tissue <- factor(data$tissue, ordered = TRUE)
p <- ggplot(data, aes(x=tissue, y=score, group=interaction(tissue, degree))) + geom_violin(width=0.7, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "darkorange2")) + geom_boxplot(width=0.08, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1 <- p + xlab("") + ylab("CpG O/E ratio") + theme(axis.title.y = element_text(family="Times",face="plain", color = "black", size = 20))
p2 <- p1 + theme(axis.text.x=element_text(family="Times",face="plain",color="black",size=rel(2.0), margin = margin(t=25))) + theme(axis.text.y=element_text(family="Times",face="plain",color="black",size=rel(2.0)))
p2 + guides(fill=F)
ggsave("CpG_OE_ratio.pdf", width = 20, height = 4)