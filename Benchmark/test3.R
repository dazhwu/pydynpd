library(plm)
library(pdynmc)
library(panelvar)

start=Sys.time()


for(i in 1:100){
  dat =read.csv('data.csv')
  pd <- pdata.frame(dat, index = c("id", "year"), drop.index = TRUE)
  z1<-pgmm(n ~ lag(n, 1:2) + w + k |lag(n, 2:4) , data=pd, effect='individual',
           model="twosteps" ,transformation='ld' , robust=TRUE)
  summary(z1, robust=TRUE)
  
  
}

print(Sys.time()-start)
summary(z1, robust=TRUE)

for (i in 1:100) {
dat =read.csv('data.csv')



m1 <- pdynmc(dat = dat, varname.i = "id", varname.t = "year",
             use.mc.diff = TRUE, use.mc.lev = FALSE, use.mc.nonlin = FALSE,
             include.y = TRUE, varname.y = "n", lagTerms.y = 2,maxLags.y=4,
             varname.reg.ex=c("w", "k"), 
             include.x=TRUE, lagTerms.reg.ex=c(0,0),
             w.mat = "iid.err", std.err = "corrected", estimation = "twostep",
             opt.meth = "none")
summary(m1)
}

print(Sys.time()-start)
# 
# 
# 
m2 <- pdynmc(dat = dat, varname.i = "id", varname.t = "year",
             use.mc.diff = TRUE, use.mc.lev = TRUE, use.mc.nonlin = FALSE,
             include.y = TRUE, varname.y = "n", lagTerms.y = 2,maxLags.y=4,
             fur.con = TRUE, fur.con.diff = TRUE, fur.con.lev = FALSE,
             varname.reg.fur = c("w", "k"), lagTerms.reg.fur = c(0,0),
             include.dum = FALSE, dum.diff = TRUE, dum.lev = FALSE, varname.dum = "year",
             w.mat = "iid.err", std.err = "corrected", estimation = "twostep",
             opt.meth = "none")
summary(m2)
mtest.fct(m2, order = 2)

start=Sys.time()
for (i in 1:100) {
dat =read.csv('data.csv')  
p1 <-pvargmm(
  dependent_vars = c("n"),
  lags = 2,
  exog_vars = c("w","k"),
  #exog_vars = c("w","k"),
  transformation = "fd",
  data = dat,
  panel_identifier = c("id", "year"),
  steps = c("twostep"),
  system_instruments = TRUE,
  max_instr_dependent_vars = 4,
  max_instr_predet_vars = 3,
  min_instr_dependent_vars = 2,
  min_instr_predet_vars = 1,
  collapse = FALSE,
  progressbar=FALSE
)

summary(p1)
}
print(Sys.time()-start)
summary(p1)


