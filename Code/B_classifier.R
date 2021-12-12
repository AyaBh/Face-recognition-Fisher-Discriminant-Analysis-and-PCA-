if (!require("OpenImageR"))
  install.packages("OpenImageR")
#require proxy to compute cosine distance
if (!require("proxy"))
  install.packages("proxy")
library("proxy")
library("OpenImageR")
#load database constants
load(file = "fisherfaces.RData")

#there's a function to test the classifiers below called testClassify

classify <- function(image) {
  # IMPORTANT: NEEDS fisherfaces.RData IN WD TO WORK
  # Predicts the person in the database that appears in the image
  # using PCA and kNN
  #
  # Args:
  #   image: character string of file in the working directory
  #
  # Returns:
  #   result: ID of the person in the database. Returns 0 if not
  #           in the database
  
  test <- t(matrix(as.vector(readImage(image))))
  
  #include eigenfaces, fda_train, train.fda, pca_train, labels
  
  
  test.scaled <- sweep(test,2,pca_train$mean)
  test.projection <- test.scaled%*%eigenfaces
  test.fda <- test.projection%*%(fda_train$P[,1:17])
  
  dmatrix <- dist(rbind(test.fda,train.fda[,1:17]),
                  method = "euclidean", diag = TRUE, upper = TRUE)
  dmatrix <- as.matrix(dmatrix)
  dmatrix <- dmatrix[1,2:(nrow(train.fda)+1)]
  sorted <- sort(dmatrix,index.return=TRUE,decreasing=FALSE)
  
  k <- 4
  labels_sel <- labels[sorted$ix][1:k]
  if (sorted$x[1] < 1598.2819119) {
    uniqv <- unique(labels_sel)
    #If they are two values with the same frecuency, it returs the first
    result <- uniqv[which.max(tabulate(match(labels_sel, uniqv)))]
  } else {
    result <- 0
  }
  
  return(result)
}

testClassify <- function(filenames, good.labels) {
  # Test the classifier on a list of images. Prints the name of the
  # file and the predicted label if incorrect.
  # Args:
  #   filenames: character vector of files in the working directory
  #
  # Returns:
  #   ccr: accuracy of the classifier on the test data
  if (length(filenames)==length(good.labels)) {
    correct <- 0
    for (i in 1:length(filenames)) {
      temp <- classify(filenames[i])
      if (temp==good.labels[i]) {
        correct <- correct + 1
      } else {
        print(paste(temp, filenames[i]))
      }
    }
    ccr <- correct/length(filenames)
    return(ccr)
  }
}