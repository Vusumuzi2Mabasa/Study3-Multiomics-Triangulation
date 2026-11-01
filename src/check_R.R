args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else '.'
library(jsonlite)
d<-fromJSON(file.path(root,'results/verification_fixture.json'),simplifyVector=FALSE)
mat<-function(a)do.call(rbind,lapply(a,unlist))
x<-mat(d$x);y<-unlist(d$y);x0<-x[,seq_len(d$p0),drop=FALSE]
Fstat<-function(y){a<-lm.fit(x0,y);b<-lm.fit(x,y);((sum(a$residuals^2)-sum(b$residuals^2))/(b$rank-a$rank))/(sum(b$residuals^2)/(length(y)-b$rank))}
fit<-lm.fit(x0,y);res<-fit$residuals/sqrt(1-rowSums(qr.Q(qr(x0))^2));weights<-mat(d$weights)
errors<-c(abs(Fstat(y)-d$F),sapply(seq_len(ncol(weights)),function(j)abs(Fstat(fit$fitted.values+res*weights[,j])-d$stars[[j]])))
z<-mat(d$z);cov<-mat(d$covariates);outcome<-unlist(d$outcome);base<-cbind(1,cov)
null<-glm.fit(base,outcome,family=binomial());logistic_errors<-numeric()
for(j in seq_len(ncol(z))){
  # Pruned MOFA scores can be 1e-130; do not rescale numerical dust into a predictor.
  if(sd(z[,j])<1e-8){p<-1}else{
    alt<-glm.fit(cbind(base,z[,j]),outcome,family=binomial())
    p<-pchisq(max(0,null$deviance-alt$deviance),1,lower.tail=FALSE)
  }
  logistic_errors<-c(logistic_errors,abs(p-d$factor_p[[j]]))
}
stopifnot(max(errors)<1e-8,max(logistic_errors)<1e-6)
write_json(list(status='passed',bootstrap_statistics=length(errors)-1,maximum_F_error=max(errors),maximum_logistic_P_error=max(logistic_errors),R=as.character(getRversion())),file.path(root,'results/R_validation.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('Independent R regression and bootstrap checks passed\n')
