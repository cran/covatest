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



#-- SYMMETRY TEST --#

sel.staz.sym <- c("DERP016", "DENW065", "DEHE051", "DETH026", "DENW063", "DENI019",
                  "DENW068", "DEHE046", "DEUB029", "DEBY047", "DETH061", "DESN049")


sp.couples.in.sym <- matrix(data = c("DERP016", "DENW065", "DEHE051", "DETH026",
                                     "DENW063", "DENI019", "DENW068", "DEHE046",
                                     "DEUB029", "DEBY047", "DETH061", "DESN049"),
                            ncol = 2, byrow = TRUE)

t.couples.in.sym <- c(1, 2)
couples.sym <- couples(sel.staz = sel.staz.sym, sp.couples.in = sp.couples.in.sym,
                       t.couples.in = t.couples.in.sym, typetest = "sym", typecode = character())

summary(couples.sym)


block.sym <- blocks(lb = 40, ls = 10, matdata = rr_13, pardata1 = 1, pardata2 = 1,
                    stpairs = couples.sym)

covabl.sym <- covablocks(stblocks = block.sym, stpairs = couples.sym,
                         typetest = "sym")

covast.sym <- covastat(matdata = rr_13, pardata1 = 1, pardata2 = 1,
                       stpairs = couples.sym, typetest = "sym")

test.sym <- covaprop(cblock = covabl.sym, cstat = covast.sym,
                     nonseptype = NULL, sign.level = 0.05)
