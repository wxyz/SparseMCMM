\name{SparseMCMM-package}
\alias{SparseMCMM-package}
\docType{package}
\title{
SparseMCMM: Estimating and testing the microbial causal mediation effect with the high-dimensional and compositional microbiome data
}
\description{
Sparse Microbial Causal Mediation Model (SparseMCMM) is designed for the high dimensional and compositional microbiome data in a treatment-microbiome-outcome causal study. SparseMCMM utilizes the linear log-contrast regression and Dirichlet regression to quantify the causal direct effect of the treatment and the causal mediation effect of the microbiome on the outcome under the counterfactual framework while addressing the compositional structure of microbiome data. Further it implements regularization techniques to handle the high-dimensional microbial mediators and identify the signature causal microbes.

SparseMCMM consists of three steps:

Step 1: Report the estimates of DE, ME and TE respectively under the sufficient causal assumptions mentioned in the Supplementary Materials, Section S1.

Step 2: Report the overall mediation test results: OME test can determine whether the overall mediation effect of microbiome is significant, and CME test can determine whether at least one individual microbe has a significant mediation effect on the outcome.

Step 3: Report the point and 95\% confidence interval estimates of $ME_j$ for each signature causal microbe identified by the regularization technique only if CME test is significant at step 2.

Consequentially, SparseMCMM provides a clear and sensible causal path analysis among treatment, compositional microbiome and outcome. Both simulation studies and real data applications showed the superb performance of SparseMCMM.}
\details{
\tabular{ll}{
Package: \tab SparseMCMM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2019-04-17\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
Chan Wang, Jiyuan Hu, Martin J. Blaser, Huilin Li.

Maintainer: Chan Wang <Chan.Wang@nyulangone.org>, Huilin Li <Huilin.Li@nyulangone.org>
}
\references{
Wang C, Hu J, Blaser M J, Li H (2019). Estimating and testing the microbial causal mediation effect with the high-dimensional and compositional microbiome data.

}
\keyword{package}
