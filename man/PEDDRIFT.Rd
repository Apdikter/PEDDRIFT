\name{PEDDRIFT}

\alias{PEDDRIFT}

\title{
PROBABILITIES OF ALLELE FREQUENCY DIFFERENCES IN DIVERGENT SELECTION LINES
}

\description{
Estimates exact probabilities of allele frequency differences between divergent selection lines.


}

\usage{
peddrift(pedigree, allele_freq=NULL, nrep=1000, rseed=0, prune=TRUE, allall=FALSE)
}

\arguments{
  \item{pedigree}{
object containing pedigree and genotype information. This file should contain the variables  tag, sire, dam, line, test (all numeric) and locus (0, 1, 2 or AA, AB, BB format).

animals with no sire or dam  (e.g., NA) are treated as founders.
}

\item{afreq}{
object (optional) containing allele frequencies in an initial founding population. If not provided, allele_freq is calculated as the (equally weighted) mean between selection lines.
}


\item{nrep}{
integer; number of replicates to be simulated under neutral (random) evolutionary drift.
}


\item{rseed}{
integer; seed used for simulations.
}

\item{prune}{
logical; if TRUE, pedigree is reduced only to genotyped animals (test = 1) and their ancestors.

(reccomended for reduced computation times).
}


\item{allall}{
logical; if TRUE, simulations wherein alleles reach fixation are excluded from calculation of test statistics.

(reccomended for rare alleles).
}


}


\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Dodds, K G and McEwan, J C (1997) Calculating exact probabilities of allele frequency differences in divergent selection lines. Proc Assoc Advmt Anim Breed Genet 12, 556-560. https://doi.org/10.13140/2.1.4568.7043
}
\author{
Zak Morrison

}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

##-- load in example data from github repo
##-- run peddrift with nrep = 1000 over 1 SNP with specified fpop allele freq 

