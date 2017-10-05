#####################################################################################################
##Function for Penalized constrained regression with glmnet                                        ##
##Double CV: inner cv for determining optimal lambda, outer for determine predictions              ##
##alpha: controls type of penalty: 0 (ridge), 0.5 (elastic net), 1 (lasso)                         ##
##nfolds: type of CV  (LOOCV: nfolds=n)                                                            ##
#####################################################################################################

glmnet.2CV<-function(X,Y,alpha=alpha,folds,nfolds,family=family,offset=NULL){

if(nrow(X)!=length(Y)) {stop("Dimensions of X and Y do not match!")}
alpha=alpha
n=nrow(X)


fit.glmnet=lapply(1:nfolds,function(i)glmnet(X[-folds[[i]],],Y[-folds[[i]]],family=family,standardize=F,alpha=alpha,offset=offset[-folds[[i]]]))
cv.fit.glmnet=lapply(1:nfolds,function(i)cv.glmnet(X[-folds[[i]],],Y[-folds[[i]]],family=family,standardize=F,alpha=alpha,offset=offset[-folds[[i]]]))

p.cv.glmnet=unlist(lapply(1:nfolds,function(i)predict(fit.glmnet[[i]],matrix(X[folds[[i]],],ncol=ncol(X)),s=cv.fit.glmnet[[i]]$lambda.min,type="response",offset=offset[folds[[i]]])))
p.cv.glmnet<-p.cv.glmnet[order(unlist(folds))]

#Calculate CV mean of the outcome#


if (family=="gaussian"){
	p0<-rep(NA,length(Y))
	glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family=family))
	for(i in 1:nfolds){
	p0[folds[[i]]]=rep(coef(glm0[[i]]))}
}
else if (family=="binomial"){
	p0<-rep(NA,length(Y))
	glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family=family))
	for(i in 1:nfolds){
	p0[folds[[i]]]=rep(exp(coef(glm0[[i]])))}
}



Q2.glmnet<-(sum((p.cv.glmnet-p0)^2)/sum((Y-p0)^2))

return(list(Q2=Q2.glmnet,p=p.cv.glmnet))
}



#####################################################################################################
##Function for asssessing Added Predictive value of X2 on top of X1                                ## 
##It requires deletion residuals (res) as input and the original Y, X1 and X2                      ##
##Double CV: inner cv for determining optimal lambda, outer for determine predictions              ##
##For performance under the null: introduce permuted event (not done internally in the function)   ##
##alpha: controls type of penalty: 0 (ridge), 0.5 (elastic net), 1 (lasso)                         ##
##nfolds: type of CV  (LOOCV: nfolds=n)                                                            ##
#####################################################################################################




AddedPredTest<-function(X1,X2,Y,res,alpha=alpha,folds,nfolds,nperm=100,family=family){



tmp<-glmnet.2CV(X2,res,alpha=alpha,folds=folds,nfolds=nfolds,family=family)
p=tmp$p
Q2=tmp$Q2


p0<-rep(NA,length(Y))
glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family=family))
for(i in 1:nfolds){
p0[folds[[i]]]=rep(coef(glm0[[i]]))}



Q2.null<-c()
for (i in 1:nperm)
	{
		
                res.null<-sample(res)
		Y.null<-(Y-res)+res.null
                p0.null<-rep(NA,length(Y.null))
                glm0<-lapply(1:nfolds,function(i)glm(Y.null[-folds[[i]]]~1,family=family))
                for(j in 1:nfolds){
                p0.null[folds[[j]]]=rep(coef(glm0[[j]]))}
		fit.null<-glmnet.2CV(X1,Y.null,alpha=alpha,folds=folds,nfolds=nfolds,family=family)
                res.null<-Y.null-fit.null$p
                fit.null<-glmnet.2CV(X2,res.null,alpha=alpha,folds=folds,nfolds=nfolds,family=family)
		Q2.null[i]<-fit.null$Q2       
	}

pvalue <- sum(Q2.null>Q2)/nperm 

return(list(Q2=Q2,pvalue=pvalue))
}





createFolds<-function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
    if (is.numeric(y)) {
        cuts <- floor(length(y)/k)
        if (cuts < 2) 
            cuts <- 2
        if (cuts > 5) 
            cuts <- 5
        y <- cut(y, unique(quantile(y, probs = seq(0, 1, length = cuts))), 
            include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(as.character(y))
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) {
            seqVector <- rep(1:k, numInClass[i]%/%k)
            if (numInClass[i]%%k > 0) 
                seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
            foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
        }
    }
    else foldVector <- seq(along = y)
    if (list) {
        out <- split(seq(along = y), foldVector)
        names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
            sep = "")
        if (returnTrain) 
            out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
    out
}


p0<-function(Y,folds,family){

	if (family=="gaussian"){
        p0<-rep(NA,length(Y))
        glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family=family))
        for(i in 1:nfolds){
        p0[folds[[i]]]=rep(coef(glm0[[i]]))}
	}
	else if (family=="binomial"){
        p0<-rep(NA,length(Y))
        glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family=family))
        for(i in 1:nfolds){
        p0[folds[[i]]]=rep(exp(coef(glm0[[i]])))}
	}
p0
}


Q2<-function(Y,p,p0){
Q2<-sum((p-p0)^2)/sum((Y-p0)^2)
Q2
}
