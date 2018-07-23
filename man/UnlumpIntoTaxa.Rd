\name{UnlumpIntoTaxa}
\alias{UnlumpIntoTaxa}
\title{
Unlump flow by clades
}
\description{
Converts the lumped description of gene flow into one line per species pair, which is what BMhyb uses in its main functions
}

\usage{
UnlumpIntoTaxa(lumped.flow)
}
\arguments{
  \item{flow}{structure of gene flow frum LumpIntoClades}
   }


\value{
	a flow table that contains the donor taxon, recipient taxon, proportion of gene flow, time from root donor, and time from root recipient.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted.
	}
