#--- load packages ---#

library("sp")
library("spacetime")
library("gstat")
library("covatest")


data(air)
ls()

if (!exists("rural")) rural = STFDF(stations, dates, data.frame(PM10 =
                                                                  as.vector(air)))
rr = rural[,"2005::2010"]

unsel = which(apply(as(rr, "xts"), 2, function(x) all(is.na(x))))

r5to10 = rr[-unsel,]

rr_13 <- r5to10[c("DEHE046","DESN049","DETH026","DENW063","DETH061","DEBY047",
                  "DENW065","DEUB029","DENW068","DENI019","DEHE051","DERP016","DENI051"),
                "2005::2006"]


#-- NON SEPARABILITY INDEX --#

C00_13<-var(rr_13[,,"PM10"]@data[[1]], na.rm = TRUE) #compute the Global Sill

#vv_13 <- variogram(PM10~1, rr_13, width=60, cutoff = 220, tlags=0:15) #estimate the spatio-temporal variogram

data(vv_13) #load the estimated spatio-temporal variogram (for more details see the help page of vv_13)

plot(vv_13, wireframe=TRUE)

nonsep.index<-sepindex(vario_st=vv_13, nt=16, ns=4, globalSill=C00_13)

boxplot(nonsep.index, ylab="Non-separability ratio")


#-- TYPE OF NON SEPARABILITY TEST --#

sel.staz.sym <- c("DERP016", "DENW065", "DEHE051", "DETH026", "DENW063", "DENI019",
                  "DENW068", "DEHE046", "DEUB029", "DEBY047", "DETH061", "DESN049")


sp.couples.in.sym <- matrix(data = c("DERP016", "DENW065", "DEHE051", "DETH026",
                                     "DENW063", "DENI019", "DENW068", "DEHE046",
                                     "DEUB029", "DEBY047", "DETH061", "DESN049"),
                            ncol = 2, byrow = TRUE)
t.couples.in.tns <- c(3, 4, 5)
couples.tns <- couples(sel.staz = sel.staz.sym, sp.couples.in = sp.couples.in.sym,
                       t.couples.in = t.couples.in.tns, typetest = "tnSep", typecode = character())

couples.tns <- setzero(x = couples.tns, zero = TRUE, value = 0)
summary(couples.tns)

block.tns <- blocks(lb=60, ls=23, matdata = rr_13, pardata1 = 1, pardata2 = 1,
                    stpairs = couples.tns)

covabl.tns <- covablocks(stblocks = block.tns, stpairs = couples.tns, typetest = "tnSep")

covast.tns <- covastat(matdata = rr_13, pardata1 = 1, pardata2 = 1, stpairs = couples.tns,
                       typetest = "tnSep")

test.tns <- covaprop(cblock = covabl.tns, cstat = covast.tns, nonseptype = 0,
                     sign.level = 0.05)
