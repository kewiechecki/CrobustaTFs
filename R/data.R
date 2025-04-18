#' PWMs for 3162 potential TF binding sites.
#'
#' A list of binding sites from ANISEED, Cis-BP, and HOMER corresponding to C. robusta transcription factors. The same TF may have multiple motifs if it is present in multiple databases.
#' 
#' @format A PWMatrixList of 3162 elements.
#' \describe{
#' 	\item{Matrix()}{The position weight matrix}
#'	\item{names()}{The unique name of each PWM.}
#'	\item{name()}{The nonunique name of each PWM.}
#'	\item{ID()}{The KH2013 Gene ID of each TF.}
#'	\item{tags()}{Additional metadata including motif source and TF family.}
#' }
#' @source \url{https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fcirobu%2FSelex_seq_Cirobu_All.zip}, \url{http://cisbp.ccbr.utoronto.ca/}, \url{http://homer.ucsd.edu/homer/download.html}
#' @seealso \link{PWMatrixList}
"CrobustaMotifs"
