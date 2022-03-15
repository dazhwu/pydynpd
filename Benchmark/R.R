library(Matrix)
library(pdynmc)
library(plm)
library(panelvar)

abdata=read.csv("data.csv")

### the following code produces an error
start=Sys.time()
for (i in 1:100){
 mc <- pdynmc(dat=abdata,varname.i = "id", varname.t = "year",
             use.mc.diff = TRUE, use.mc.lev = TRUE, use.mc.nonlin = FALSE,
             include.y = TRUE, varname.y = "n", lagTerms.y = 2, include.x = TRUE,
             varname.reg.pre = c("w", "k"), lagTerms.reg.pre = c(0,0), maxLags.reg.pre = c(3,3),
             include.dum = FALSE, dum.diff = FALSE, dum.lev = FALSE,
             w.mat = "iid.err", std.err = "corrected", estimation = "twostep",
             opt.meth = "none"
)
summary(mc)
}
print(Sys.time()-start)

### the following code produces results inconsistent with other packages
start=Sys.time()
for (i in 1:100){
 mc <- pdynmc(dat=abdata,varname.i = "id", varname.t = "year",
             use.mc.diff = TRUE, use.mc.lev = TRUE, use.mc.nonlin = FALSE,
             include.y = TRUE, varname.y = "n", lagTerms.y = 2, include.x = TRUE,
             varname.reg.pre = c("w", "k"), lagTerms.reg.pre = c(0,0), maxLags.reg.pre = c(3,3),
             include.dum = FALSE, dum.diff = FALSE, dum.lev = FALSE,
             w.mat = "iid.err", std.err = "corrected", estimation = "twostep",
             opt.meth = "none"
)
summary(mc)
}
print(Sys.time()-start)


### panelvar
start=Sys.time()
for (i in 1:100){
ex3_abdata <-pvargmm(
  dependent_vars = c("n"),
  lags = 2,
  predet_vars = c("w"),
  exog_vars=c("k"),
  transformation = "fd",
  data = abdata,
  panel_identifier = c("id", "year"),
  steps = c("twostep"),
  system_instruments = TRUE,
  max_instr_dependent_vars = 3,
  max_instr_predet_vars = 3,
  min_instr_dependent_vars = 1L,
  min_instr_predet_vars = 1L,
  collapse = FALSE
)
summary(ex3_abdata)
}
print(Sys.time()-start)


######  plm  ################

start=Sys.time()
for (i in 1:100){

pd <- pdata.frame(abdata, index = c("id", "year"), drop.index = TRUE)
z1<-pgmm(n ~ 1+ lag(n, 1:2) + w + k |lag(n, 2:4) + lag(w, 1:3), data=pd, effect='individual',
         model="twosteps" ,transformation='ld' , fsm='FULL')
summary(z1, robust=TRUE)

}
print(Sys.time()-start)
