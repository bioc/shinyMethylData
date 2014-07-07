\name{summary.tcga.raw}
\alias{summary.tcga.raw}
\docType{data}
\title{
Example dataset of input data for shinyMethyl
}
\description{
Extracted data from 369 TCGA Head and Neck Cancer DNA methylation samples. The extracted data serve as an example dataset for the package shinyMethyl. Original samples are from 450k methylation arrays, and were obtained from The Cancer Genome Atlas (TCGA). 310 samples are from tumor, 50 are matched normals and 9 are technical replicates of a control cell line. }
\usage{data(summary.tcga.raw)}
\format{
A \code{list} containing the necessary information to launch a shinyMethyl session. See the links below for more details on the data.
}
\references{
The Cancer Genome Atltas (TCGA) Head and Neck Cancer dataset:
http://cancergenome.nih.gov/cancersselected/headandneck
}
\seealso{
See \code{\link[shinyMethyl]{shinySummarize}} for details on how to perform the data extraction. 
See \code{\link[shinyMethyl]{runShinyMethyl}} for how to launch a shinyMethyl session.
 }
\examples{
data(summary.tcga.raw)	
\dontrun{
runShinyMethyl(summary.tcga.raw)
}
}
\keyword{datasets}

