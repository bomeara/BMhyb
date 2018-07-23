\name{AdjustFlow}
\alias{AdjustFlow}
\title{
	AdjustFlow
	}

\description{
	Randomly deletes taxa of hybrid origin or hybridization events to make the network consistent with our restrictions (that there is only one hybridization event leading to any taxon). Some important things to note: 1) This is changing history. You say there's a hybridization event, and it instead deletes it. That's not good. Or it deletes taxa, also not good. [but, it's the best available in this program -- and after all, you were happy using trees with NO hybridization events for years...]. 2) It does it randomly: so you run an analysis once after doing AdjustFlow, and get one answer, but if you do it again, you might get a different answer.
	}

\usage{
	AdjustFlow(data,phy,flow, remove=c("taxa", "events"))
	}

\arguments{
\item{data}{continuous trait data containing species information in vector format.}
\item{phy}{a tree in phylo class.}
\item{flow}{a structure of gene flow.}
\item{remove}{randomly delete taxa of hybrid origin ("taxa") or hybridization events ("events")}
	}

\details{
By choosing "taxa" the final phylogeny and data vector will be small, so those are modified. By choosing "events" those are untouched. In both cases, the flow data.frame is modified.
}

\value{
	A list containing final data, final phylogeny, final flow, and the number of attempts required (equivalent to how many taxa / events were deleted).
 	}

\author{
	Dwueng-Chwuan Jhwueng & Brian O'Meara
	}
