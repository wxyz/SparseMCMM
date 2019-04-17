

SparseMCMM=function(Treatment,otu.com,outcome,covariates=NULL,covariate.fix=NULL,dirichlet.penalty=seq(0,1,0.1),
                            lm.penalty1=seq(0,1,0.1),lm.penalty2=seq(0,2,0.2),
                            low.bound1=NULL,up.bound1=NULL,low.bound2=NULL,up.bound2=NULL,
                            max.iter=3000,num.per=NULL,bootstrap=NULL)
{


  ########### parameter estimations

  alpha.estimation=alpha.estimates(Treatment,otu.com,outcome,covariates,
                           penalty.lambda1=lm.penalty1,penalty.lambda2=lm.penalty2,
                           low.bound=low.bound1,up.bound=up.bound1)


  beta.estimation=beta.estimates(Treatment,otu.com,covariates,penalty.lambda=dirichlet.penalty,
                                 low.bound=low.bound2,up.bound=up.bound2)



  otu.num=ncol(otu.com)

  CausalEffect=CausalE(otu.com,alpha.estimation,beta.estimation,covariate.fix)

  Effect.estimates=CausalEffect



  if(!is.null(num.per)) {

    ME=abs(Effect.estimates[[1]][2])
    CME= sum(Effect.estimates[[2]]^2)

    sample.num=length(Treatment)
    Z=log(otu.com)
    X=cbind(rep(1,sample.num),Treatment,covariates,Z,Z*Treatment)

    outcome.predict=c(X%*%alpha.estimation)

    residual=outcome-outcome.predict


    for(rr in 1:num.per)
    {
      outcome.per=outcome.predict+sample(residual,length(residual))

      Treatment.per=sample(Treatment,length(Treatment))



      alpha.estimation.per=alpha.estimates(Treatment,otu.com,outcome.per,covariates,
                                       penalty.lambda1=lm.penalty1,penalty.lambda2=lm.penalty2,
                                       low.bound=low.bound1,up.bound=up.bound1)

      beta.estimation.per=beta.estimates(Treatment,otu.com,covariates,penalty.lambda=dirichlet.penalty,
                                     low.bound=low.bound2,up.bound=up.bound2)



      CausalEffect1=CausalE(otu.com,alpha.estimation,beta.estimation.per,covariate.fix)
      CausalEffect2=CausalE(otu.com,alpha.estimation.per,beta.estimation,covariate.fix)
      CausalEffect3=CausalE(otu.com,alpha.estimation.per,beta.estimation.per,covariate.fix)



      ME=c(ME,max(abs(CausalEffect1[[1]][2]),abs(CausalEffect2[[1]][2]),abs(CausalEffect3[[1]][2])))
      CME= c(CME,max(sum(CausalEffect1[[2]]^2),sum(CausalEffect2[[2]]^2),sum(CausalEffect3[[2]]^2)))
    }



    rr=cbind(ME,CME)

    p.values=apply(rr,2,function(x){sum(x>=x[1])/length(x)})

    names(p.values)=c("OME","CME")

  }


  if(!is.null(bootstrap)) {

    individual.ME=NULL

    for(bb in 1:bootstrap)
    {

      in.b=c(sample(which(Treatment==0),length(which(Treatment==0)),replace = TRUE),
                   sample(which(Treatment==1),length(which(Treatment==1)),replace = TRUE))

      Treatment.per=Treatment[in.b]
      outcome.per=outcome[in.b]
      otu.com.per=otu.com[in.b,]
      if(!is.null(covariates)) covariates.per=covariates[in.b,] else covariates.per=NULL

      alpha.estimation=alpha.estimates(Treatment.per,otu.com.per,outcome.per,covariates.per,
                                       penalty.lambda1=lm.penalty1,penalty.lambda2=lm.penalty2,
                                       low.bound=low.bound1,up.bound=up.bound1)

      beta.estimation=beta.estimates(Treatment.per,otu.com.per,covariates.per,penalty.lambda=dirichlet.penalty,
                                     low.bound=low.bound2,up.bound=up.bound2)

      CausalEffect=CausalE(otu.com.per,alpha.estimation,beta.estimation,covariate.fix)

      individual.ME=rbind(individual.ME,CausalEffect[[2]])

    }


    me.95upper=Effect.estimates[[2]]+sqrt(diag(cov(individual.ME))/bootstrap)*1.96
    me.95lower=Effect.estimates[[2]]-sqrt(diag(cov(individual.ME))/bootstrap)*1.96

    me.ind=rbind(Effect.estimates[[2]],me.95upper,me.95lower)

    rownames(me.ind)=c("Compontent-wiseME","95%UpperCI","95%LowerCI")

    if(is.null(colnames(otu.com))) colnames(me.ind)=1:otu.num else colnames(me.ind)=colnames(otu.com)

    me.ind=me.ind[,(sign(me.95lower)*sign(me.95upper))==1]



  }


if(!is.null(num.per)) Effect.estimates$Test= p.values
if(!is.null(bootstrap)) Effect.estimates$'Compontent-wise ME 95% CI'= me.ind


return(Effect.estimates)

}




