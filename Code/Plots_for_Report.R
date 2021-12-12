library(ggplot2)
pca <- pca(data)
data.scaled <- sweep(data,2,pca$mean)
class.num    <- length(unique(labels)) - 1
eigenfaces <- t(data.scaled)%*%pca$P[,1:class.num]
data.projection <- data.scaled%*%eigenfaces
data.pca.lab <- data.frame(data.projection, labels)
fda <- fda(data.pca.lab)
data.fda <- data.frame(data.projection%*%fda$P)
short <- c(1,2,3,4,5,6,7,8,9,10,12,18,19)
data.fda$hair <- rep(NA,150)

for (i in 1:150) {
  if (labels[i] %in% short) {
    data.fda$hair[i] <- "short"
  } else {
    data.fda$hair[i] <- "long"
  }
}

data.new = as.data.frame(data.fda)
data.new$labels = factor(labels)
data.new$names = factor(names)

ggplot(data.new, aes(LDA1, LDA2,color = labels, fill = labels))+
  geom_point(alpha=0.5, position="identity", bins = 50)+
  theme(legend.position = "bottom")

ggplot(data.new, aes(data.new[,1], color = labels, fill = labels))+
  geom_histogram(alpha=0.5, position="identity", bins = 50)+
  theme(legend.position = "bottom")+
  facet_wrap(vars(labels), ncol = 5)

ggplot(data.new, aes(LDA1, LDA2, color=labels, label=labels))+
  geom_label()+
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

tres<-ggplot(data.new, aes(LDA2, LDA3, color=hair, label=labels))+
  geom_label()+
  theme(legend.position = "right", axis.text = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

ggplot(data.new[labels%in%c(4,14,24),], aes(LDA1, LDA2, color=hair, label=names))+
  geom_label()+
  theme(legend.position = "right", axis.text = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

colors <- c("#999999", "#E69F00")
colors <- colors[as.numeric(iris$Species)]
scatterplot3d::scatterplot3d(data.new[,1:3], color = colors)

library(plotly)
fig <- plot_ly(data.new, x = ~LDA1, y = ~LDA2, z = ~LDA3, color = ~hair, colors = c('#BF382A', '#0C4B8E'))
