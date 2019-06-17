## ---- include=FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(serac)

## ----eval=FALSE----------------------------------------------------------
#  getwd()

## ----eval=FALSE----------------------------------------------------------
#  setwd("/myProject/Age depth model")

## ----eval=FALSE----------------------------------------------------------
#  dir.create(file.path(getwd(), 'Cores'), showWarnings = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  dir.create(file.path(paste(getwd(),'/Cores',sep=""), 'serac_example_ALO09P12'), showWarnings = FALSE)
#  

## ----eval=FALSE----------------------------------------------------------
#  ?serac_example_ALO09P12
#  write.table(x = serac_example_ALO09P12, file = paste(getwd(),'/Cores/serac_example_ALO09P12/serac_example_ALO09P12.txt',sep=""),col.names = T, row.names = F,sep="\t")
#  
#  # Including proxy data for this core too
#  write.table(x = serac_example_ALO09P12_proxy, file = paste(getwd(),'/Cores/serac_example_ALO09P12/serac_example_ALO09P12_proxy.txt',sep=""),col.names = T, row.names = F,sep="\t")
#  

## ----eval=FALSE----------------------------------------------------------
#  serac(name="serac_example_ALO09P12",coring_yr=2009)

## ----eval=FALSE----------------------------------------------------------
#  serac(name="serac_example_ALO09P12",coring_yr=2009,model=c("CFCS"),plotphoto=FALSE,minphoto=c(0),maxphoto=c(210),plot_Pb=T,plot_Am=T,plot_Cs=T,Cher=c(30,40),Hemisphere=c("NH"),NWT=c(51,61),sedchange=c(75.5),plot_Pb_inst_deposit=T,inst_deposit=c(20,28,100,107,135,142,158,186),suppdescriptor=TRUE,descriptor_lab=c("Ca/Fe"),historic_d=c(20,28,100,107,135,142,158,186),historic_a=c(1994,1920,1886,1868),historic_n=c("sept1 994 flood","1920 flood","1886 flood","1868 flood ?"), min_yr=c(1750),dmax=c(180), plotpdf=TRUE,preview=F)
#  

## ----eval=FALSE----------------------------------------------------------
#  mymodel1 <- serac(name="serac_example_ALO09P12",coring_yr=2009)
#  
#  # See the objects created
#  names(mymodel1)
#  
#  # Visualize the input data
#  mymodel1$data
#  
#  # Save the plot with custom parameters, e.g., Courier font and grey scale without editing default colors.
#  pdf(paste0(getwd(),"/Cores/serac_example_ALO09P12/serac_example_ALO09P12_custom.pdf"),width = 10, height = 5, family="Courier", colormodel = "grey")
#  mymodel1$plot
#  dev.off()
#  

