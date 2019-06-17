## ---- include=FALSE---------------------------------------------------------------------
stopifnot(require(knitr))
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L,width = 90)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = FALSE,
  dev = "png",
  dpi = 150,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "90%",
  fig.align = "center"
)

library(devtools)
devtools::install_github("rosalieb/serac")
library(serac)


## ----eval=FALSE-------------------------------------------------------------------------
#  getwd()

## ----eval=FALSE-------------------------------------------------------------------------
#  setwd("/myProject/Age depth model")

## ----eval=TRUE--------------------------------------------------------------------------
dir.create(file.path(getwd(), 'Cores'), showWarnings = FALSE)

## ----eval=TRUE--------------------------------------------------------------------------
dir.create(file.path(paste(getwd(),'/Cores',sep=""), 'serac_example_ALO09P12'), showWarnings = FALSE)


## ----eval=TRUE--------------------------------------------------------------------------
?serac_example_ALO09P12
write.table(x = serac_example_ALO09P12, file = paste(getwd(),'/Cores/serac_example_ALO09P12/serac_example_ALO09P12.txt',sep=""),col.names = T, row.names = F,sep="\t")

# Including proxy data for this core too
write.table(x = serac_example_ALO09P12_proxy, file = paste(getwd(),'/Cores/serac_example_ALO09P12/serac_example_ALO09P12_proxy.txt',sep=""),col.names = T, row.names = F,sep="\t")


