\name{core_metadata}
\alias{core_metadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Core metadata
}
\description{
Question the user on core International Geo Sample Number, sampling number, etc.
}
\usage{
core_metadata(name = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
Name of the core, given using quotes. Use preferably the published name of the core for traceability.
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
%% ~put references to the literature/web site here ~
}
\author{
Rosalie Bruel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (name = c())
{
    answer = NA
    answer2 = NA
    answer3 = NA
    rm(answer, answer2, answer3)
    message("You're about to enter some informations on the method used to measure radioelements on your core.\nType Y to proceed, N to abord (and turn the argument archive_metadata to FALSE if you just \nran the serac function) so this message doesn't show up again.")
    answer <- readline()
    if (answer == "n")
        answer = "N"
    if (answer == "y")
        answer = "Y"
    while (TRUE) {
        if (answer == "N") {
            message("\nTurn archive_metadata to FALSE and re run your code line.\n")
            break
        }
        else {
            answer = "N"
            if (length(list.files(paste(getwd(), "/Cores/", name,
                "/", sep = ""), pattern = "serac_metadata_suppmetadata*",
                full.names = TRUE)) == 1) {
                message("\nA file with metadata already exists in this folder. Do you want to:\n    - replace it: Tape 1 in the console\n    - edit it: Tape 2 in the console")
                answer2 <- readline()
            }
            if (exists("answer2") && answer2 == "2")
                suppmetadata <- read.delim(list.files(paste(getwd(),
                  "/Cores/", name, "/", sep = ""), pattern = "serac_metadata_suppmetadata*",
                  full.names = TRUE), header = T, sep = "")
            b = 1
            cat("\nInternational Geo Sample Number (IGSN) in the System for Earth Sample Registration Database.\nUnique ID of your core, see www.geosamples.org \n  (e.g., EDYSAN004)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                ISGN <- paste(suppmetadata[b, 2])
            if (answer3 == "Y") {
                message("  Enter ISGN:")
                ISGN <- readline()
            }
            b = b + 1
            cat("\nWhen was the core sampled? (YYYY-MM-DD) \n  (e.g., 2015-11-14)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                sample_date <- paste(suppmetadata[b, 2])
            if (answer3 == "Y") {
                message("  Enter sample date (YYYY-MM-DD):")
                sample_date <- readline()
            }
            b = b + 1
            cat("\nWhat is the latitude (y) of the coring location in WGS84? \n  (e.g., 46.020976)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                coring_coordinates_y <- paste(suppmetadata[b,
                  2])
            if (answer3 == "Y") {
                message("  Enter latitude:")
                coring_coordinates_y <- readline()
            }
            b = b + 1
            cat("\nWhat is the longitude (x) of the coring location in WGS84? \n  (e.g., 6.768192)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                coring_coordinates_x <- paste(suppmetadata[b,
                  2])
            if (answer3 == "Y") {
                message("  Enter longitude:")
                coring_coordinates_x <- readline()
            }
            b = b + 1
            cat("\nWhat was the coring method? \n  (e.g., gravity corer, piston corer, percussion, etc.)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                coring_method <- paste(suppmetadata[b, 2])
            if (answer3 == "Y") {
                message("  Enter coring method:")
                coring_method <- readline()
            }
            b = b + 1
            cat("\nWhat was the subsampling method at the laboratory? \n  (e.g., calibrated volumetric sampler, extrusion, etc.)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                laboratory_subsampling_method <- paste(suppmetadata[b,
                  2])
            if (answer3 == "Y") {
                message("  Enter subsampling method:")
                laboratory_subsampling_method <- readline()
            }
            b = b + 1
            cat("\nIn which lab measurements where carried out? \n  (e.g., Edytem, FR)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                laboratory <- paste(suppmetadata[b, 2])
            if (answer3 == "Y") {
                message("  Enter laboratory name:")
                laboratory <- readline()
            }
            b = b + 1
            cat("\nInstrument type: \n  (e.g., alpha spectrometry, well-type germanium detector, P-type germanium detector)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                instrument_type <- paste(suppmetadata[b, 2])
            if (answer3 == "Y") {
                message("  Enter instrument type:")
                instrument_type <- readline()
            }
            b = b + 1
            cat("\nStart of the measurements (YYYY-MM-DD):\n  (e.g., 2018-07-01)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                measurement_startdate <- paste(suppmetadata[b,
                  2])
            if (answer3 == "Y") {
                message("  Enter start of the measurements:")
                measurement_startdate <- readline()
            }
            b = b + 1
            cat("\nEnd of the measurements (YYYY-MM-DD):\n  (e.g., 2018-07-02)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                measurement_enddate <- paste(suppmetadata[b,
                  2])
            if (answer3 == "Y") {
                message("  Enter end of the measurements:")
                measurement_enddate <- readline()
            }
            b = b + 1
            cat("\nAny additional comments?\n  (e.g., 210Pb background reached, first centimetres lost during coring, etc.)")
            if (exists("answer2") && answer2 == 2) {
                message(paste("  Current information in your file is: ",
                  suppmetadata[b, 2], ". Replace? (Y/N)", sep = ""))
                answer3 <- readline()
                if (answer3 == "n")
                  answer3 = "N"
                if (answer3 == "y")
                  answer3 = "Y"
            }
            else {
                answer3 = "Y"
            }
            if (answer3 == "N")
                additional_comments <- paste(suppmetadata[b,
                  2])
            if (answer3 == "Y") {
                message("  Enter additional comments:")
                additional_comments <- readline()
            }
            suppmetadata <- matrix(c("Parameter", "Answer given",
                "ISGN", ISGN, "sample_date", sample_date, "coring_coordinates_y",
                coring_coordinates_y, "coring_coordinates_x",
                coring_coordinates_x, "coring_method", coring_method,
                "laboratory_subsampling_method", laboratory_subsampling_method,
                "measurement_laboratory", laboratory, "instrument_type",
                instrument_type, "measurement_startdate", measurement_startdate,
                "measurement_enddate", measurement_enddate, "additional_comments",
                additional_comments), ncol = 2, byrow = T)
            cat("_______________________________________________________")
            message(paste("You're all set! Please contact the authors if you would \nlike to see any other metadata included.\n",
                sep = ""))
            suppmetadata2 <- as.data.frame(suppmetadata[-1, 2])
            rownames(suppmetadata2) <- suppmetadata[-1, 1]
            colnames(suppmetadata2) <- suppmetadata[1, 2]
            print(suppmetadata2)
            cat("_______________________________________________________")
            write.table(x = suppmetadata, file = paste(getwd(),
                "/Cores/", name, "/serac_metadata_suppmetadata.txt",
                sep = ""), col.names = F, row.names = F)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
