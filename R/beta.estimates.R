#require(Compositional)
#require(nloptr)

beta.estimates=function(Treatment,otu.com,covariates=NULL,penalty.lambda=seq(0,1,0.1),
                        low.bound=NULL,up.bound=NULL,max.iter=3000)
{
  sample.num=nrow(otu.com);p=ncol(otu.com)

  if(is.null(covariates)) q=0 else q=ncol(covariates)


  beta0.intial=diri.est(otu.com[Treatment==0,],type = "mle")$param

  beta.initial=c(beta0.intial,rep(0,(1+q)*p))


  BIC=estimates=NULL

  for(lambda.dirichlet in penalty.lambda)
  {
    qq=nloptr( x0=beta.initial,
               toOpt2,
               lb=low.bound,
               ub=up.bound,
               opts =list( "algorithm" = "NLOPT_LD_SLSQP",
                           "xtol_rel" = 1.0e-4,
                           "maxeval"= max.iter),
               otu.com=otu.com,Treatment=Treatment,
               covariates=covariates,lambda.dirichlet=lambda.dirichlet)

    if(qq$status %in% c(1,3,4)){

      estimates1=round(qq$solution,3)

      estimates=rbind(estimates,estimates1)

      BIC=c(BIC,log(sample.num)*sum(estimates1!=0)-
              2*likelihood2.org(estimates1,Treatment=Treatment,covariates=covariates,otu.com=otu.com))
    }}

  if(length(BIC)==0) beta.estimates=beta.initial else beta.estimates=estimates[which.min(BIC),]

  return(beta.estimates)

}
