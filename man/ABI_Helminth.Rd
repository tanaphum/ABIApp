% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABI_Helminth.R
\name{ABI_Helminth}
\alias{ABI_Helminth}
\title{A tool for helminth species delimitation at various taxonomic levels}
\usage{
ABI_Helminth(
  distance = 0,
  Fastafile = NULL,
  fastaSelect1 = 1,
  fastaSelect2 = 2,
  group = "NAS",
  marker = "18S rRNA"
)
}
\arguments{
\item{distance}{Number >= 0 (From our database the number sh ould be 0.000 - 0.909)}

\item{Fastafile}{directory of Fasta file}

\item{fastaSelect1}{Selection of 1st taxon from Fasta file (number or text(label))}

\item{fastaSelect2}{Selection of 2nd taxon from Fasta file (number or text(label))
("18S rRNA","28S rRNA","ITS1","ITS2","COI","COII","cytB","NAD1","12S rRNA","16S rRNA")}

\item{group}{group of Helminth ("NAS","NS","NT","TR","TRD","CE") /n
"NAS" is "Nematode (Ascaridida and Spirurida)" /n
"NS"  is "Nematode (Strongylida)") /n
"NT"  is "Nematode (Trichocephalida)" /n
"TR"  is "Trematode (Plagiorchiida)" /n
"TRD" is "Trematode (Diplostomida)" /n
"CE"  is "Cestode"}

\item{marker}{Helminth Genetic Markers}
}
\value{
Plot of ggplot
}
\description{
A tool for helminth species delimitation at various taxonomic levels
}
\examples{
ABI_Helminth()
ABI_Helminth(0.06)
ABI_Helminth(0.02,"NS","18S rRNA")
ABI_Helminth(distance = 0.5,group = "CE",marker = "ITS2")

##### Fasta file #####
# Select with numbers
ABI_Helminth(Fastafile="dir/fastaFile.fasta","NT","18S rRNA"
fastaSelect1 = 1, fastaSelect2 = 2
)
# Select with number and text
ABI_Helminth(Fastafile="dir/fastaFile.fasta","NT","18S rRNA"
fastaSelect1 = "Label 1", fastaSelect2 = 2
)
# Select with texts
ABI_Helminth(Fastafile="dir/fastaFile.fasta","NT","18S rRNA"
fastaSelect1 = "Label 1", fastaSelect2 = "Label 2"
)

}
