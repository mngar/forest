#' Function to cross validate genomic selection by rrblup
#' @param x genotype matrix
#' @param pheno phenotype vector
#' @param cv number of cross validation
#' @param projName name of the project ("character")
rrblup <- function(x, pheno, cv, part, projName){
#generacion de index
train <- round(part/100*nrow(x), digits = 0)          # part == porcentaje pob entrenamiento
test <- nrow(x)-train
if(isTRUE(test+train == nrow(x)) == TRUE) {
xtrain <- rep(1,train)
xtest <- rep(2,test)
xcv <- c(xtrain, xtest)
index <- matrix(nrow = nrow(x), ncol = cv, NA)
for (i in 1:cv) {
index[,i] <- sample(xcv)
}
} else {
print("train+test es distinto de nrow(x), adecuar valores")
}
################################################################################
datas <- cbind(x, pheno)
sets1  <- index
datas2 <- population
datas3 <- feno
ntest <- test
for(fold in 1:cv){
 correlations <- matrix(NA,1,6)
 datos        <- matrix(NA,1,6)
################################################################################
  itrain <- which(sets1[,fold]==1)
  itest  <- which(sets1[,fold]==2)
  test   <- datas[itest,]
  train  <- datas[itrain,]
  Xtest  <- test[,-ncol(test)]
  Ytest  <- test[,ncol(test)]
  Xtrain <- train[,-ncol(test)]
  Ytrain <- train[,ncol(test)]
################################################################################
indice <- 1
for (columna1 in 1:(ncol(Xtrain)-1)){
    for (columna2 in (columna1+1):ncol(Xtrain)){
        dd <- sum(abs(Xtrain [,columna1]- Xtrain [,columna2]))
    if(dd==0){
        indice <- c(indice,columna1)
}}}
##  RR-BLUP       ##############################################################
  X1     <- rbind(Xtrain,Xtest)
  X1     <- X1[,-indice]
  y1     <- c(Ytrain,Ytest)
  yNa1   <- y1
  train1 <- 1:nrow(Xtrain)
  f      <-  nrow(Xtrain)+1
  pred1  <- f:nrow(X1)
  ans    <- mixed.solve(y=y1[train1],Z=Xtrain) #By default K = I
  intercepto <- rep(ans$beta,ntest)
###
efectoRR  <- ans$u
efectoRR2 <- abs(efectoRR)
num <- length(efectoRR)
rr  <- cbind(fold,1:num,efectoRR,efectoRR2)
newdata <- rr[order(-efectoRR2),]
yTestHat   <- intercepto + Xtest %*% ans$u
accuracyRR <- cor(Ytest,yTestHat)
prediccion <- paste(fold,Ytest,yTestHat, sep = ",")
write.table(prediccion, file = paste0(projName, "_observado_vs_predicho.csv"),row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
write.table(newdata,file= paste0(projName, "_SALIDA_EFECTOS_RR.csv"),row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
write.table(accuracyRR,file= paste0(projName, "_salida_PRECISION_SG.csv"),row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
    }
write.table(index,file= paste0(projName, "_Index.csv"))
}




