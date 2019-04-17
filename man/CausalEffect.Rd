\name{CausalEffect}
\alias{CausalEffect}
\title{
Calculating causal DE, ME and component-wise ME estimates based on models (1)-(2)
}
\description{
This function calculates the causal DE of treatment, ME and component-wise ME of microbiome on the outcome under the counterfactual framework with the estimated coefficients in models (1)-(2).
}

\usage{

CausalE(otu.com,alpha.estimation,beta.estimation,covariate.fix)
}


\arguments{
  \item{otu.com}{
  A n*p numeric matrix containing compositional microbiome data. Each row represents a subject, and each column represents a taxon (given the rank, for example, the genus rank) or an OTU. The row sum equals 1.
}

  \item{alpha.estimation}{
  The numeric vector containing estimated coefficients from alpha.estimates function
}

  \item{beta.estimation}{
  The numeric vector containing estimated coefficients from beta.estimates function
}

  \item{covariate.fix}{
  An optional vector containing the given values for covariates. Default=NULL.
}

}


\value{
A list which contains two elements:
  \describe{
  \item{Esitmated Causal Effects}{It has three values: direct effect (DE), mediation effect (ME) and total effect (TE)}
  \item{Estimated component-wise Mediation Effects}{A numeric vector containing component-wise MEs for all mediators}
  }
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
# require(dirmult)
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
# beta.estimation=beta.estimates(Treatment,otu.com,covariates)
#
# alpha.estimation=alpha.estimates(Treatment,otu.com,outcome,covariates)
#
# CausalE(otu.com,alpha.estimation,beta.estimation,covariate.fix=c(0,0))

}
\keyword{CausalEffect}

