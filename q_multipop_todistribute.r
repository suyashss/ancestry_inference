eps=1e-4

gt_texttonum<-function(gtstr){
	x=strsplit(gtstr,":")[[1]][1]
	ans=0
	if(x=="0/1" || x=="0|1" || x=="1/0" || x=="1|0"){
		ans=1
	}
	else if(x=="1|1" || x=="1/1"){
		ans=2
	}
	else{
		ans=0
	}
	ans
}

boundf<-function(f){
	f1=f
	f1[f1<eps]=eps;f1[f1>1-eps]=1-eps
	f1
}

ll<-function(x,g,fmat){ #fmat=JxK
	qvec=x#c(x,1-sum(x))
	qmat=matrix(qvec,nrow=length(qvec),ncol=1,byrow=F) #Kx1
	fbar=boundf(c(fmat%*%qmat))
	-1*sum(g*log(fbar)+(2-g)*log(1-fbar))
}

gradll<-function(x,g,fmat){ #fmat=JxK
	qvec=x#c(x,1-sum(x))
	qmat=matrix(qvec,nrow=length(qvec),ncol=1,byrow=F)# Kx1
	fbar=boundf(c(fmat%*%qmat)) #Jx1
	gmat=matrix(g,ncol=length(qvec),nrow=length(g),byrow=F) # JxK
	fbarmat=matrix(fbar,nrow=length(fbar),ncol=length(qvec),byrow=F)
	gradvec=c(colSums( fmat*gmat/fbarmat+(2-gmat)*(1-fmat)/(1-fbarmat)))
	-1*gradvec
}

issumq1<-function(x,g,fmat){
	sum(x)-1
}

issumq1_g<-function(x,g,fmat){
	rep(1,length(x))
}

init_em<-function(g,fmat){
	K=dim(fmat)[2]
	xold=rep(1/K,K)
	xnew=xold
	gmat=matrix(g,ncol=K,nrow=length(g),byrow=F) # JxK
	for(i in 1:200 ){
		qmat=matrix(xold,nrow=K,ncol=1,byrow=F)# Kx1
		fbar=boundf(c(fmat%*%qmat)) #Jx1
		fbarmat=matrix(fbar,nrow=length(fbar),ncol=K,byrow=F)
		qmat2=matrix(xold,nrow=dim(fmat)[1],ncol=K,byrow=T)
		amat=qmat2*fmat/fbarmat
		bmat=qmat2*(1-fmat)/(1-fbarmat)
		xnew=colSums(gmat*amat+(2-gmat)*bmat)/(2*length(fbar))
		xnew=xnew/sum(xnew)
		#cat("EM:",i,xold,"--->",xnew,"=",sum(xnew),"\n")
		xold=boundf(xnew)
	}
	#cat("Initial guess is",xnew,"\n")
	boundf(xnew)
}

findq<-function(g,fmat){
	K=dim(fmat)[2]
	local_opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1e-7)
	opts=list("algorithm"="NLOPT_LD_SLSQP",print_level=0,maxeval=10000,xtol_rel=1e-7,ftol_abs=1e-10)
	x0=init_em(g,fmat)
	solobj=nloptr(x0=x0,eval_f=ll,eval_grad_f=gradll,eval_g_eq=issumq1,eval_jac_g_eq=issumq1_g,lb=rep(eps,K),ub=rep(1-eps,K),g=g,fmat=fmat,opts=opts)
	bestq=solobj$solution
	#print(solobj)
	#print(bestq)
	bestq
}

args<-commandArgs(T)
if(length(args)<2){
	stop("Usage:Rscript q_multipop_todistribute.r frequencyfile vcffile\n")
}

library("plyr")
#library("data.table")
library("nloptr")
#library("ade4")

frqfile<-args[[1]]
vcffile<-args[[2]]

frqtable<-read.table(frqfile,header=T)
vcffull<-read.table(vcffile,header=F,stringsAsFactors=F)
colnames(vcffull)<-c('CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','IND')

vcf=vcffull[,c(1,2,4,5,10)]

ob=join(frqtable,vcf,type="left")
#head(ob)
#print(dim(frqtable))
#print(dim(ob))

poplist=c("AFR","ASN","AMR","SAN","EUR")
freq_subset=1-as.matrix(ob[,poplist])
gtvec_text=ob$IND
gtvec_text[is.na(gtvec_text)]="0|0"

gtvec=sapply(gtvec_text,gt_texttonum)
#head(gtvec)

qvec=findq(gtvec,freq_subset)
names(qvec)<-poplist
print(qvec)
