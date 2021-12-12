setwd("C:/Users/Alex/Downloads/Training/Training")
if (!require("OpenImageR"))
  install.packages("OpenImageR")
if (!require("proxy"))
  install.packages("proxy")
if (!require("tictoc"))
  install.packages("tictoc")

library("OpenImageR")
library("proxy")
library("tictoc")
# short(1,2,3,4,5,6,7,8,9,10,12,19)
names <- list.files(pattern = "jpg")

Img1 = readImage(names[1])
data <- matrix(0, length(names),prod(dim(Img1)))
dim(data)
labels <- vector()

for (i in 1:length(names)){
  data[i,] <- as.vector(readImage(names[i]))
  if (nchar(names[i]) > 7) {
    labels[i] <- as.integer(substr(names[i],1,2))
  } else {
    labels[i] <- as.integer(substr(names[i],1,1))
  }
}

pca <- function(data) {
  # Computes useful features to use PCA to reduce dimensionality
  #
  # Args:
  #   data: data whose sample covariance is to be calculated.
  #
  # Returns:
  #   mean: Vector containing the mean of the columns of the data
  #   P: Matrix containing the eigenvectors of the dispersion matrix of the data 
  #   D: vector containing the variance explained by each principal component
  
  mean <- colMeans(data)
  data.scaled <- sweep(data, 2, mean)
  pc.num <- min(c(nrow(data), ncol(data)))
  if (pc.num==ncol(data)) {
    Sigma <- t(data.scaled)%*%data.scaled/(nrow(data.scaled)-1)
  } else {
    Sigma <- data.scaled%*%t(data.scaled)/(nrow(data.scaled)-1)
  }
  Eigen <- eigen(Sigma)
  Eigenvalues <- Eigen$values
  Prop.Var <- Eigenvalues/sum(Eigenvalues)
  Eigenvectors <- Eigen$vectors
  
  
  pc.colnames <- vector()
  for (i in 1:pc.num) {
    pc.colnames <- c(pc.colnames, paste0("PC",i))
  }
  colnames(Eigenvectors) <- pc.colnames
  newList <- list("mean" = mean, "P" = Eigenvectors, "D" = Prop.Var)
  return(newList)
}

fda <- function(data) {
  # Computes useful features to use Fisher Discriminant Analysis
  #
  # Args:
  #   data: data to be projected, has a labels column on its last column
  #
  # Returns:
  #   mean: Vector containing the mean of the columns of the data
  #   P: Matrix containing the eigenvectors of the appropriate matrix 
  #   D: vector containing the variance explained by each Fisher discriminant dimension
  
  #separate the elements from our data
  labels <- data[,ncol(data)]
  data <- data[,-ncol(data)]
  data$labels=factor(labels)
  class.num <- length(unique(data$labels)) - 1
  
  #let's calculate the mean 
  Means = sapply(levels(data$labels), 
                 FUN = function(x){
                   colMeans(data[data$labels==x, 1:class.num])
                 })
  M=colMeans(data[,1:class.num])
  
  #lets calculate S.W
  sw=0
  for (i in 1:class.num) {
    sw=sw+cov(data[data$labels==i, 1:class.num])*(table(data$labels)[1]-1)
  }
  
  #let's calculate S.B
  sb=0
  for (i in 1:class.num) {
    sb=sb+(table(data$labels)[i])*(Means[,i]-M)%*%t(Means[,i]-M)
  }
  
  eig = eigen(solve(sw)%*%sb)
  eigenvectors = eig$vectors
  eigenvalues = eig$values
  Prop.Var <- eigenvalues/sum(eigenvalues)
  
  ld.colnames <- vector()
  
  for (i in 1:class.num) {
    ld.colnames <- c(ld.colnames, paste0("LDA",i))
  }
  colnames(eigenvectors) <- ld.colnames
  newList <- list("mean" = M, "P" = eigenvectors, "D" = Prop.Var)
  return(newList)
}

# data.new = as.data.frame(fisherfaces[,1])
# data.new$labels = data$labels
# ggplot(data.new, aes(data.new[,1],color = labels, fill = labels))+
#   geom_histogram(alpha=0.5, position="identity", bins = 50)+
#   theme(legend.position = "bottom")
# 
# ggplot(data.new, aes(data.new[,1], color = labels, fill = labels))+
#   geom_histogram(alpha=0.5, position="identity", bins = 50)+
#   theme(legend.position = "bottom")+
#   facet_wrap(vars(labels), ncol = 5)

findK <- function(obs, train.data, method, maxK) {
  # Returns which label the knn algorithm would give to a certain
  # observation
  #
  # Args:
  #   obs: observation whose distances will be computed.
  #   train.data: data that whose distances to the observation will
  #               be computed.
  #   method: distance metric character string.
  #
  # Returns:
  #   result: Vector containing the label whose k (column number)
  #           are nearest of our observation. The last column is the
  #           real label of the observation.
  
  label <- obs[,1]
  train.label <- train.data[,1]
  dmatrix <- dist(rbind(obs[2:length(obs)],train.data[,2:ncol(train.data)]),
                  method = method, diag = TRUE, upper = TRUE)
  dmatrix <- as.matrix(dmatrix)
  dmatrix <- dmatrix[1,2:(nrow(train.data)+1)]
  sorted <- sort(dmatrix,index.return=TRUE,decreasing=FALSE)
  result <- rep(0, maxK+2)
  result[maxK+1] <- label
  for (k in 1:maxK) {
    labels_sel <- train.label[sorted$ix][1:k]
    if (sorted$x[1] < 1598.2819119) { 
      uniqv <- unique(labels_sel)
      #If they are two values with the same frecuency, it returs the first
      result[k] <- uniqv[which.max(tabulate(match(labels_sel, uniqv)))]
    } else {
      result[k] <- 0
    }
  }
  longest_correct_distance_index <- which(train.label[sorted$ix] %in% label)
  result[maxK+2] <- sorted$x[longest_correct_distance_index[5]]
  return(result)
}

testMethod <- function(test.data, train.data, method, maxK) {
  # Computes how accurate is the knn method for a certain data and
  # distance method.
  #
  # Args:
  #   test.data: data where observations to compute the distance
  #              will be taken.
  #   train.data: data that whose distances to the observation will
  #               be computed.
  #   method: distance metric character string.
  #
  # Returns:
  #   result: vector containing the mean accuracy for every k on 
  #           all the observations in the test set. Last element
  #           is the mean of all k's.
  
  m <- findK(test.data[1,], train.data, method, maxK)
  for (i in 2:nrow(test.data)) {
    m <- rbind(m, findK(test.data[i,], train.data, method, maxK))
  }
  meanKs <- numeric()
  for (i in 1:(ncol(m)-2)) {
    meanKs <- c(meanKs, sum(m[,i]==m[,ncol(m)-1])/nrow(m))
  }
  m <- rbind(m, c(meanKs,0,0))
  m[nrow(m),ncol(m)-1] <- mean(m[nrow(m),1:(ncol(m)-2)])
  m[nrow(m),ncol(m)] <- max(m[1:(nrow(m)-1),ncol(m)])
  return(m[nrow(m),])
}

plotGenerator <- function(th, maxK) {
  # Generates a bar plot that describes the accuracy of the model
  # using a certain threshold
  #
  # Args:
  #   th: number of eigenvectors retained 
  #
  # Returns:
  #   max_correct_distance: the maximum distance from the farthest correct observation
  #   p: ggplot plot
  
  #add columns to identify each observation
  data.wLabel <- data.frame(labels = labels, x=data)
  data.wLabel <- data.wLabel[order(labels),]
  data.wLabel <- data.frame(index=as.integer(rownames(data.wLabel)), data.wLabel)
  
  #matrices to store the repeated cross-validation
  euclideanMatrix <- vector()
  manhattanMatrix <- vector()
  cosineMatrix <- vector()
  
  #Cross-validation
  for (i in 1:10) {
    #Splitting the data
    library(splitstackshape)
    #k-fold using just the index and the label
    a <- stratified(data.wLabel[,1:2],"labels",1,bothSets =TRUE)
    b <- stratified(a$SAMP2,"labels",1,bothSets =TRUE)
    c <- stratified(b$SAMP2,"labels",1,bothSets =TRUE)
    d <- stratified(c$SAMP2,"labels",1,bothSets =TRUE)
    e <- stratified(d$SAMP2,"labels",1,bothSets =TRUE)
    one <- a$SAMP1[,1]
    two <- b$SAMP1[,1]
    three <- c$SAMP1[,1]
    four <- d$SAMP1[,1]
    five <- e$SAMP1[,1]
    six <- e$SAMP2[,1]
    kfold_idx <- list(one,two,three,four,five,six)
    
    data.noIndex <- data.wLabel[,-1] #remove indexes
    
    #matrix to store accuracy for a repetition
    euclidean <- vector()
    manhattan <- vector()
    cosine <- vector()
    
    #k-fold
    for (i in 1:6) {
      train_data <- as.matrix(data.noIndex[-unlist(kfold_idx[i]),-1])
      test_data <- as.matrix(data.noIndex[unlist(kfold_idx[i]),-1])
      train_labels <- data.noIndex[-unlist(kfold_idx[i]),1]
      test_labels <- data.noIndex[unlist(kfold_idx[i]),1]
      
      pca_train <- pca(train_data)
      train.scaled <- sweep(train_data,2,pca_train$mean)
      class.num    <- length(test_labels) - 1
      
      eigenfaces <- t(train.scaled)%*%pca_train$P[,1:class.num]
      
      train.projection <- train.scaled%*%eigenfaces
      train.pca.lab <- data.frame(train.projection, train_labels)
      
      fda <- fda(train.pca.lab)
      
      train.fda <- train.projection%*%(fda$P[,1:th])
      train.new <- data.frame(labels = train_labels,x = train.fda)
      
      test.scaled <- sweep(test_data,2,pca_train$mean)
      test.projection <- test.scaled%*%eigenfaces
      test.fda <- test.projection%*%(fda$P[,1:th])

      test.new <- data.frame(labels = test_labels,x = test.fda)
      
      euclidean <- rbind(euclidean, testMethod(test.new, train.new, method = "euclidean", maxK))
      manhattan <- rbind(manhattan, testMethod(test.new, train.new, method = "manhattan", maxK))
      cosine <- rbind(cosine, testMethod(test.new, train.new, method = "cosine", maxK))
    }
    euclideanMatrix <- rbind(euclideanMatrix, colMeans(euclidean))
    manhattanMatrix <- rbind(manhattanMatrix, colMeans(manhattan))
    cosineMatrix <- rbind(cosineMatrix, colMeans(cosine))
  }
  
  #resultsMatrix SPECIFY NUMBER OF DISTANCES AND NAMES
  results <- rbind(colMeans(euclideanMatrix), colMeans(manhattanMatrix), colMeans(cosineMatrix))
  max_correct_distance <- results[,ncol(results)]
  rownames(results) <- c("Euclidean", "Manhattan", "Cosine")
  results <- as.data.frame(results[,-c(ncol(results), ncol(results)-1)])
  colnames(results) <- 1:maxK
  
  #to define matrix to plot
  k <- maxK
  distances <- 3
  distanceNameVector <- character()
  for (i in 1:distances) {
    distanceNameVector <- c(distanceNameVector, rep(rownames(results)[i],k))
  }
  #matrix for using with ggplot
  resultsPlot <- data.frame(Distance=distanceNameVector,
                            K=rep(1:k, distances),
                            Value=rep(0,k*distances))
  l <- 1
  for (i in 1:distances) {
    for (j in 1:k) {
      resultsPlot[l,3] <- round(results[i,j],3)
      l <- l + 1
    }
  }
  
  resultsPlot <- as.data.frame(resultsPlot)
  colnames(resultsPlot) <- c("Distance", "K", "Value")
  resultsPlot$Value <- as.numeric(resultsPlot$Value)
  
  # To get the percentage of variance for the plot
  pca <- pca(data)
  data.scaled <- sweep(data,2,pca$mean)
  class.num    <- length(unique(labels)) - 1
  eigenfaces <- t(data.scaled)%*%pca$P[,1:class.num]
  data.projection <- data.scaled%*%eigenfaces
  data.pca.lab <- data.frame(data.projection, labels)
  fda <- fda(data.pca.lab)
  
  library(ggplot2)
  p <- ggplot(resultsPlot)+
    aes(fill=Distance, y=Value, x=K) +
    ylim(0,1)+
    ylab("Accuracy")+
    labs(title=paste("Retaining ", th, " eigenvectors, ", round(sum(fda$D[1:th])*100,1),
                     "% of the variance",sep = ""), subtitle = "Out-of-sample accuracy")+
    geom_bar(position="dodge", stat="identity")+
    geom_text(aes(label=Value), position=position_dodge(width=0.9), vjust=-0.25)+
    scale_x_continuous(breaks = 1:maxK)
  newList <- list("max_correct_distance"=max_correct_distance, "p"=p)
  return(newList)
}

# To see how many eigenvectors we need for a specific variance
# min(which(cumsum(fda$D)>0.95))
tic()
nine <- plotGenerator(17,10)
toc()
nine$p
nine$max_correct_distance
