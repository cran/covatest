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



#-- TEST ON THE TYPE OF MODEL --#


#In couples.mod the spatio-temporal lags to be used for the tests on the type of model
# are provided. This object will be used for the tests on the Product-Sum,
# Integrated Product and Gneiting class of models

sel.staz.mod <- c("DERP016", "DENW065", "DENW063", "DEHE046", "DEUB029",
                  "DETH061", "DENW068", "DETH026", "DENI051")

sp.couples.in.mod <- matrix(data = c("DERP016", "DENW065", "DENW063", "DEHE046",
                                     "DEUB029", "DETH061", "DEHE046", "DENW063",
                                     "DERP016", "DENW068", "DETH026", "DENI051",
                                     "DEUB029", "DETH061", "DENI051", "DETH061",
                                     "DERP016", "DEUB029"),
                            ncol = 2, byrow = TRUE)

t.couples.in.mod <- c(1, 2, 3)

couples.mod <- couples(sel.staz = sel.staz.mod, sp.couples.in =
                         sp.couples.in.mod, t.couples.in = t.couples.in.mod,
                       typetest = "productSum", typecode = character())

zero.index <- matrix(data=c(3, 7, 6, 7, 9, 7), ncol=2, byrow = TRUE)

couples.mod <- setzero(x = couples.mod, zero = FALSE, index = zero.index, value = 0)


#In block.mod, the blocks to be used for the estimation of the covariance matrix
# for the tests on the type of model are provided.
#This object will be used for the tests on the Product-Sum, Integrated Product
# and Gneiting class of models
block.mod <- blocks(lb = 60, ls = 10, matdata = rr_13, pardata1 = 1,
                    pardata2 = 1, stpairs = couples.mod)



#### Product-Sum Model###
covabl.ps <- covablocks(stblocks = block.mod, stpairs = couples.mod,
                        typetest = "productSum")

covast.ps <- covastatM(matdata = rr_13, pardata1 = 1, pardata2 = 1,
                       stpairs = couples.mod, typetest = "productSum", beta.data = NULL)

test.ps <- covaprop(cblock = covabl.ps, cstat = covast.ps, nonseptype = NULL,
                    sign.level = 0.05)


#### Integrated Product Model ###
#For this test the spatio-temporal lags used for the product-sum model (couples.mod)
# and the same blocks (block.mod) have been considered

covabl.ip <- covablocks(stblocks = block.mod, stpairs = couples.mod,
                        typetest = "intProduct")

covast.ip <- covastatM(matdata = rr_13, pardata1 = 1, pardata2 = 1,
                       stpairs = couples.mod, typetest = "intProduct", beta.data = NULL)

test.ip <- covaprop(cblock = covabl.ip, cstat = covast.ip, nonseptype = NULL,
                    sign.level = 0.05)


#### Gneiting Model###
#For this test the spatio-temporal lags used for the product-sum model (couples.mod)
# and the same blocks (block.mod) have been considered

covabl.gn <- covablocks(stblocks = block.mod, stpairs = couples.mod,
                        typetest = "gneiting")

covast.gn <- covastatM(matdata = rr_13, pardata1 = 1, pardata2 = 1,
                       stpairs = couples.mod, typetest = "gneiting", beta.data = 1)

test.gn <- covaprop(cblock = covabl.gn, cstat = covast.gn, nonseptype = NULL,
                    sign.level = 0.05)
