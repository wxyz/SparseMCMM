#######################Parameter estimation for the Dirichlet regression model

######The log likelihood function for the Dirichlet regression

likelihood2.org=function(beta,Treatment,covariates,otu.com){

  num.otu=ncol(otu.com); sample.num=nrow(otu.com);

  beta0=beta[1:num.otu];

  if(is.null(covariates)) {q=0;betaX=NULL} else {q=ncol(covariates);

  betaX=beta[(1+num.otu):((q+1)*num.otu)];
  betaX=matrix(betaX,nrow=num.otu,ncol=q)}


  betaT=beta[-(1:((q+1)*num.otu))];

  ##### betaC is  num.otu*q

  X=cbind(rep(1,sample.num),covariates,Treatment) #n*(1+q+p)

  b=cbind(beta0,betaX,betaT) #num.otu*(1+q+p)

  g=exp(X %*% t(b))	# n * num.otu

  ### log likelihood

  sum(lgamma(rowSums(g)) + rowSums((g-1) * log(otu.com) - lgamma(g)))
}


######The log likelihood function for the Dirichlet regression with lasso penalty for variable selection

likelihood2=function(beta,Treatment,covariates,otu.com,lambda.dirichlet){

  num.otu=ncol(otu.com);

  beta.pen=beta[-(1:num.otu)]

  ### -log likelihood+ penalty

  -likelihood2.org(beta,Treatment,covariates,otu.com)+lambda.dirichlet*(sum(abs(beta.pen)))
}


######The gradient function for the log likelihood with lasso penalty

gradient2 = function(beta,Treatment,covariates,otu.com,lambda.dirichlet) {

  num.otu=ncol(otu.com); sample.num=nrow(otu.com);

  beta0=beta[1:num.otu];

  if(is.null(covariates)) {q=0;betaX=NULL} else {q=ncol(covariates);

  betaX=beta[(1+num.otu):((q+1)*num.otu)];
  betaX=matrix(betaX,nrow=num.otu,ncol=q)}


  betaT=beta[-(1:((q+1)*num.otu))];

  X=cbind(rep(1,sample.num),covariates,Treatment) #n*(1+q+p)

  b=cbind(beta0,betaX,betaT) #num.otu*(1+q+p)

  g=exp(X %*% t(b))	# n * num.otu

  S=t((digamma(rowSums(g))  - digamma(g) + log(otu.com)) * g) %*% X  ## 1+q+p

  beta.pen=c(betaX,betaT)

  der=rep(0,length(beta.pen))

  if(sum(beta.pen!=0)>0) der[beta.pen!=0]=lambda.dirichlet*beta.pen[beta.pen!=0]/abs(beta.pen[beta.pen!=0])

  c(-S)+c(rep(0,num.otu),der)

}

#### function to optimize - a list of objective and gradient
toOpt2 = function(beta,Treatment,covariates,otu.com,lambda.dirichlet){
  list(objective=likelihood2(beta,Treatment,covariates,otu.com,lambda.dirichlet),
       gradient=gradient2(beta,Treatment,covariates,otu.com,lambda.dirichlet))
}
