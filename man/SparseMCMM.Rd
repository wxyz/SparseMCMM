\name{SparseMCMM}
\alias{SparseMCMM}
\title{A main function in SparseMCMM framework}
\description{

This function gives estimated DE, estimated ME, estimated TE, and estimated component-wise MEs, statistical significances of OME and CME with permutation procedure, as well as, 95\% CIs of component-wise MEs for the causal mediators with bootstrapping procedure, based on models (1)-(2).
}


\usage{
SparseMCMM(Treatment,otu.com,outcome,covariates=NULL,covariate.fix=NULL,
           dirichlet.penalty=seq(0,1,0.1),lm.penalty1=seq(0,1,0.1),lm.penalty2=seq(0,2,0.2),
          low.bound1=NULL,up.bound1=NULL,low.bound2=NULL,up.bound2=NULL,
           max.iter=3000,num.per=NULL,bootstrap=NULL)
}



\arguments{
  \item{Treatment}{
A numeric vector of the binary treatment (takes the value 1 if it is assigned to the treatment group and takes the value 0 if assigned to the control group) with length = sample size (n).
}
  \item{otu.com}{
  A n*p numeric matrix containing compositional microbiome data. Each row represents a subject, and each column represents a taxon (given the rank, for example, the genus rank) or an OTU. The row sum equals 1.
}

  \item{outcome}{
A numeric vector of the continuous outcome with length = sample size (n).
}


  \item{covariates}{
  An optional matrix containing covariates which need to be adjusted in the model.Each row represents a subject, and each column represents a covariate. Default=NULL.
}

  \item{covariate.fix}{
  An optional vector containing the given values for covariates. Default=NULL.
}

  \item{dirichlet.penalty}{
A numeric vector constaing the candidated tunning parameters used in beta.estimates function .Default=seq(0,1,0.1).
}

  \item{lm.penalty1}{
A numeric vector constaing the candidated tunning parameters for the first penalty function in the alpha.estimates.Default=seq(0,1,0.1).
}

  \item{lm.penalty2}{
A numeric vector constaing the candidated tunning parameters for the second penalty function in the alpha.estimates.Default=seq(0,2,0.2).
}

  \item{low.bound1}{
A numeric vector with lower bounds of the controls used in alpha.estiates function, Default=NULL, there is no lower bound.
}

  \item{up.bound1}{
A numeric vector with upper bounds of the controls used in alpha.estiates function, Default=NULL, there is no upper bound.
}

  \item{low.bound2}{
A numeric vector with lower bounds of the controls used in beta.estiates function, Default=NULL, there is no lower bound.
}

  \item{up.bound2}{
A numeric vector with upper bounds of the controls used in beta.estiates function, Default=NULL, there is no upper bound.
}

  \item{max.iter}{
An integer value, the maximum number of iterations, Default=3000.
}

  \item{num.per}{
An integer value, the number of permuations. statistical significances of tests TME and CME are calculated based on these permutations. Default=NULL, No calculation for hypothesis test.
}


  \item{bootstrap}{
An integer value, the number of bootstrappings. 95\% CIs of component-wise MEs are calculated based on these bootstrappings. Default=NULL, no calculation for 95\% CIs.
}

}

\value{

A list which contains four elements:
  \describe{
  \item{Esitmated Causal Effects}{It has three values: direct effect (DE), mediation effect (ME) and total effect (TE)}
  \item{Estimated component-wise Mediation Effects}{A numeric vector containing component-wise MEs for all mediators}

  \item{Test}{P-values for tests OME and CME if num.per is not NULL, otherwise, no this item.}

  \item{Compontent-wise ME 95\% CI}{A numeric matrix if bootstrap is not NULL, otherwise, no this item.}

  }}

\references{
Wang C, Hu J, Blaser M J, Li H (2019). Estimating and testing the microbial causal mediation effect with the high-dimensional and compositional microbiome data.
}
\author{
Chan Wang, Jiyuan Hu, Martin J. Blaser, Huilin Li.
}

\examples{

# ########### generation data
# ### Sample size and number of mediators
# sample.num=100
# otu.num=10
#
# ###Treatment
# Treatment=rep(c(0,1),each=sample.num/2)
# ######Two covariates
# covariates=cbind(sample(c(1,0),sample.num,replace = TRUE),rnorm(sample.num))
#
# ### parameters
# beta0=c(0.6, -0.3, 0.8, -1.4, -1.2, -1.4, -1.3, -1.0, -0.2, 0.6)
# betaT=rep(0,otu.num)
# betaT[c(1,3)]=c(0.4,0.2)
# betaX=matrix(0,otu.num,2)
#
# alpha0=0
# alphaT=1
# alphaZ=alphaC=rep(0,otu.num)
# alphaZ[c(1,3)]=c(0.7,-0.7)
# alphaC[c(1,3)]=c(0.15,-0.15)
# alphaX=c(0,0)
#
#
# ############Microbiome data
# library(dirmult)
# X=cbind(rep(1,sample.num),covariates,Treatment) #n*(1+q+p)
# b=cbind(beta0,betaX,betaT) #num.otu*(1+q+p)
# gamma.simu=exp(X %*% t(b))	# n * num.otu
# otu.com=t(apply(gamma.simu,1,rdirichlet,n=1))
#
# ##################Outcome  data
# X=cbind(rep(1,sample.num),Treatment,covariates,log(otu.com),log(otu.com)*Treatment)
# b=c(alpha0,alphaT,alphaX,alphaZ,alphaC)
# outcome=c(b%*%t(X)+rnorm(sample.num,mean = 0, sd =1))
#
#
#
# ##### SparseMCMM function
# SparseMCMM(Treatment,otu.com,outcome,covariates,covariate.fix=c(0,0),
#           num.per=10,bootstrap=10)

}
\keyword{SparseMCMM}







