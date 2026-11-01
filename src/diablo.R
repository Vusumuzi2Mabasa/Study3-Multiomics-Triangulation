# Official mixOmics supervised signatures, trained independently of MOFA/eQTL.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else '.'
suppressPackageStartupMessages(library(mixOmics));library(jsonlite)
pilot<-'--pilot' %in% args
ba<-function(y,p)mean(c(mean(p[y==0]==0),mean(p[y==1]==1)))
readmat<-function(path)as.matrix(read.csv(path,header=FALSE,check.names=FALSE))
prepare<-function(a,b,ca,cb){
  xa<-cbind(1,ca[,3:6,drop=FALSE]);xb<-cbind(1,cb[,3:6,drop=FALSE]);oa<-ob<-list()
  for(name in names(a)){
    coefficient<-lm.fit(xa,a[[name]])$coefficients;coefficient[!is.finite(coefficient)]<-0
    aa<-a[[name]]-xa%*%coefficient;bb<-b[[name]]-xb%*%coefficient
    good<-apply(aa,2,sd)>1e-7
    oa[[name]]<-aa[,good,drop=FALSE];ob[[name]]<-bb[,good,drop=FALSE]
  }
  oa$clinical<-ca;ob$clinical<-cb
  colnames(oa$clinical)<-colnames(ob$clinical)<-paste0('C_',1:ncol(ca))
  list(train=oa,test=ob)
}
fitmodel<-function(x,y,k,keep){
  design<-matrix(.1,length(x),length(x));diag(design)<-0
  keepx<-lapply(x,function(a)rep(min(keep,ncol(a)),k));keepx$clinical<-rep(ncol(x$clinical),k)
  withCallingHandlers(block.splsda(x,factor(y,levels=c(0,1)),ncomp=k,keepX=keepx,design=design,scale=TRUE,max.iter=1000),warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
}
predictclass<-function(model,x,k)as.integer(as.character(predict(model,newdata=x)$WeightedVote$centroids.dist[,k]))
selected<-function(model,block){w<-model$loadings[[block]];as.integer(sub(paste0('^',block,'_'),'',rownames(w)[rowSums(abs(w))>0]))-1L}
files<-list.files(file.path(root,'results/cohorts'),pattern='[.]json$',full.names=TRUE)
if(pilot)files<-files[grepl('(null|convergent)_240_00[.]json$',files)]
for(file in files){
  id<-sub('[.]json$','',basename(file));out<-file.path(root,'results/diablo',paste0(id,'.json'))
  if(file.exists(out))next
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE);meta<-fromJSON(file);set.seed(meta$seed+181)
  warnings_seen<-character();began<-proc.time()[3]
  result<-tryCatch({
    folder<-file.path(root,'data/cohorts',id)
    blocks<-function(label){x<-lapply(c('genotype','RNA'),function(name){v<-readmat(file.path(folder,paste0(label,'_',name,'.csv')));colnames(v)<-paste0(name,'_',seq_len(ncol(v)));v});names(x)<-c('genotype','RNA');x}
    a<-blocks('train');b<-blocks('test');ca<-readmat(file.path(folder,'train_covariates.csv'));cb<-readmat(file.path(folder,'test_covariates.csv'))
    ya<-scan(file.path(folder,'train_y.csv'),quiet=TRUE);yb<-scan(file.path(folder,'test_y.csv'),quiet=TRUE)
    folds<-integer(length(ya));for(group in c(0,1)){ix<-sample(which(ya==group));folds[ix]<-rep(1:3,length.out=length(ix))}
    candidates<-list()
    for(k in c(1,2))for(keep in c(3,6)){
      scores<-numeric(3)
      for(f in 1:3){tr<-folds!=f;te<-!tr;p<-prepare(lapply(a,function(x)x[tr,,drop=FALSE]),lapply(a,function(x)x[te,,drop=FALSE]),ca[tr,,drop=FALSE],ca[te,,drop=FALSE]);fit<-fitmodel(p$train,ya[tr],k,keep);scores[f]<-ba(ya[te],predictclass(fit,p$test,k))}
      candidates[[length(candidates)+1]]<-list(k=k,keep=keep,BA=mean(scores))
    }
    best<-candidates[[which.max(sapply(candidates,function(x)x$BA))]];p<-prepare(a,b,ca,cb)
    model<-fitmodel(p$train,ya,best$k,best$keep);prediction<-predictclass(model,p$test,best$k)
    stopifnot(all(prediction %in% c(0,1)))
    stability<-list()
    if(meta$replicate==0)for(boot in 1:5){
      ix<-unlist(lapply(c(0,1),function(group){positions<-which(ya==group);sample(positions,length(positions),replace=TRUE)}))
      pp<-prepare(lapply(a,function(x)x[ix,,drop=FALSE]),b,ca[ix,,drop=FALSE],cb)
      bm<-fitmodel(pp$train,ya[ix],best$k,best$keep)
      stability[[boot]]<-list(genes=selected(bm,'RNA'),variants=selected(bm,'genotype'))
    }
    list(status='ok',scenario=meta$scenario,n=meta$n,replicate=meta$replicate,seed=meta$seed,k=best$k,keep=best$keep,candidates=candidates,genes=selected(model,'RNA'),variants=selected(model,'genotype'),test_class=prediction,balanced_accuracy=ba(yb,prediction),stability=stability,warnings=unique(warnings_seen),converged=!any(grepl('converg',warnings_seen,ignore.case=TRUE)),R=as.character(getRversion()),mixOmics=as.character(packageVersion('mixOmics')),seconds=unname(proc.time()[3]-began))
  },error=function(e)list(status='failed',scenario=meta$scenario,n=meta$n,replicate=meta$replicate,error=conditionMessage(e)))
  write_json(result,out,auto_unbox=TRUE,pretty=TRUE,digits=16);cat(id,result$status,'\n')
}
