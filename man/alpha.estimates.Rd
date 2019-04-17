\name{alpha.estimates}
\alias{alpha.estimates}
\title{
Parameter estimation for the linear log-contrast model
}
\description{
This function gives the numerical estimates of the parameters in model (1) by minimizing the penalized sum of squared residuals  measuring the discrepancy between the observed and predicted outcome with sequential quadratic programming (SQP) method. With the heredity constraint, two penalties are incorporated to preform variable selection in high dimensional linear regression
with interaction (optimization problem (7)). The tunning parameter is determined by Bayesian information criterion (BIC).
}
\usage{

alpha.estimates(Treatment,otu.com,outcome,covariates=NULL,
                         penalty.lambda1=seq(0,1,0.1),penalty.lambda2=seq(0,2,0.2),
                         low.bound=NULL,up.bound=NULL,max.iter=3000)
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

  \item{penalty.lambda1}{
A numeric vector constaing the candidated tunning parameters for the first penalty function.Default=seq(0,1,0.1).
}

  \item{penalty.lambda2}{
A numeric vector constaing the candidated tunning parameters for the second penalty function.Default=seq(0,2,0.2).
}


  \item{low.bound}{
A numeric vector with lower bounds of the controls, Default=NULL, there is no lower bound.
}

  \item{up.bound}{
A numeric vector with upper bounds of the controls, Default=NULL, there is no upper bound.
}

  \item{max.iter}{
An integer value, the maximum number of iterations, Default=3000.
}

}


\value{
A numeric vector, the parameters'estimates. There is no optimal numerical estimates if return 0.
}

\references{
Wang C, Hu J, Blaser M J, Li H (2019). Estimating and testing the microbial causal mediation effect with the high-dimensional and compositional microbiome data.
}
\author{
Chan Wang, Jiyuan Hu, Martin J. Blaser, Huilin Li.
}

\examples{
# require(SparseMCMM)
#
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
# ########
# alpha.estimates(Treatment,otu.com,outcome,covariates)

}
\keyword{alpha.estimates}

