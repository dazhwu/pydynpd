
The objective of the package is similar to the following open-source packages: <br>
Package | Language
--- | --- 
plm | R
panelvar | R
pdynmc | R

To compare pydynpd with similar packages, we performed several benchmark tests. More specifically, for each package we run 100 times to estimate the same model with the same data. For verification, the tests also include Stata package xtabond2. This is because xtabond2 is the most popular package in dynamic panel model.

Please note for several models, we do not include pdynmc because pdynmc kept on reporting error messages and we couldn't figure out how to make it work.

> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Deepin 20.4

Matrix products: default
BLAS/LAPACK: /opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so


numpy.show_config()
blas_armpl_info:
  NOT AVAILABLE
blas_mkl_info:
    libraries = ['mkl_rt', 'pthread', 'mkl_rt']
    library_dirs = ['/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64']
    define_macros = [('SCIPY_MKL_H', None), ('HAVE_CBLAS', None)]
    include_dirs = ['/opt/intel/compilers_and_libraries/linux/mkl/include']
blas_opt_info:
    libraries = ['mkl_rt', 'pthread', 'mkl_rt']
    library_dirs = ['/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64']
    define_macros = [('SCIPY_MKL_H', None), ('HAVE_CBLAS', None)]
    include_dirs = ['/opt/intel/compilers_and_libraries/linux/mkl/include']

## Test 1 Difference GMM

plm:

library(plm)
dat =read.csv('data.csv')
pd <- pdata.frame(dat, index = c("id", "year"), drop.index = TRUE)
z1<-pgmm(n ~ lag(n, 1:2) + w + k |lag(n, 2:4) , data=pd, effect='individual',
model="twosteps" ,transformation='d' , robust=TRUE)
           
summary(z1)

Unbalanced Panel: n = 140, T = 7-9, N = 1031

Number of Observations Used: 611
Residuals:
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-0.671621 -0.035738  0.000000 -0.004808  0.037571  0.560347 

Coefficients:
              Estimate Std. Error z-value  Pr(>|z|)    
lag(n, 1:2)1  0.391272   0.204642  1.9120 0.0558786 .  
lag(n, 1:2)2 -0.089571   0.071560 -1.2517 0.2106824    
w            -0.498584   0.129817 -3.8407 0.0001227 ***
k             0.435816   0.072260  6.0312 1.628e-09 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Sargan test: chisq(15) = 35.41221 (p-value = 0.0021474)
Autocorrelation test (1): normal = -1.469485 (p-value = 0.1417)
Autocorrelation test (2): normal = -0.02836377 (p-value = 0.97737)
Wald test for coefficients: chisq(4) = 313.7221 (p-value = < 2.22e-16)

pdynmc
m2 <- pdynmc(dat = dat, varname.i = "id", varname.t = "year",
             use.mc.diff = TRUE, use.mc.lev =FALSE, use.mc.nonlin = FALSE,
             include.y = TRUE, varname.y = "n", lagTerms.y = 2,maxLags.y=4,
             fur.con = TRUE, fur.con.diff = TRUE, fur.con.lev = FALSE,
             varname.reg.fur = c("w", "k"), lagTerms.reg.fur = c(0,0),
             w.mat = "iid.err", std.err = "corrected", estimation = "twostep",
             opt.meth = "none")
summary(m2)
mtest.fct(m2, order = 2)
Dynamic linear panel estimation (twostep)
Estimation steps: 2

Coefficients:
     Estimate Std.Err.rob z-value.rob Pr(>|z.rob|)    
L1.n  0.39127     0.20464       1.912      0.05588 .  
L2.n -0.08957     0.07156      -1.252      0.21057    
L0.w -0.49858     0.12982      -3.841      0.00012 ***
L0.k  0.43582     0.07226       6.031      < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

 19 total instruments are employed to estimate 4 parameters
 17 linear (DIF) 
 2 further controls (DIF) 
 no time dummies 
 
J-Test (overid restrictions):  35.41 with 15 DF, pvalue: 0.0021
F-Statistic (slope coeff):  313.72 with 4 DF, pvalue: <0.001
F-Statistic (time dummies):  no time dummies included in estimation
> mtest.fct(m2, order = 2)

	Arellano and Bond (1991) serial correlation test of degree 2

data:  2step GMM Estimation
normal = -0.034135, p-value = 0.9728
alternative hypothesis: serial correlation of order 2 in the error terms

panelvar
p1 <-pvargmm(
  dependent_vars = c("n"),
  lags = 2,
  exog_vars = c("w","k"),
  #exog_vars = c("w","k"),
  transformation = "fd",
  data = dat,
  panel_identifier = c("id", "year"),
  steps = c("twostep"),
  system_instruments = FALSE,
  max_instr_dependent_vars = 4,
  max_instr_predet_vars = 3,
  min_instr_dependent_vars = 2,
  min_instr_predet_vars = 1,
  collapse = FALSE,
  progressbar=FALSE
)

summary(p1)

---------------------------------------------------
Dynamic Panel VAR estimation, two-step GMM 
---------------------------------------------------
Transformation: First-differences 
Group variable: id 
Time variable: year 
Number of observations = 611 
Number of groups = 140 
Obs per group: min = 4 
               avg = 4.364286 
               max = 6 
Number of instruments = 19 

===================
        n          
-------------------
lag1_n   0.3913    
        (0.2046)   
lag2_n  -0.0896    
        (0.0716)   
w       -0.4986 ***
        (0.1298)   
k        0.4358 ***
        (0.0723)   
===================
*** p < 0.001; ** p < 0.01; * p < 0.05

---------------------------------------------------
Instruments for  equation
 Standard
  FD.(w k)
 GMM-type
  Dependent vars: L(2, 4)
  Collapse =  FALSE 
---------------------------------------------------

Hansen test of overid. restrictions: chi2(15) = 35.41 Prob > chi2 = 0.002
(Robust, but weakened by many instruments.)

df = pd.read_csv("data.csv")
command_str = 'n L(1:2).n w k  | gmm(n, 2:4)  iv(w k) | nolevel'
mydpd = regression.abond(command_str, df, ['id', 'year'])
    
Dynamic panel-data estimation, two-step difference GMM
 Group variable: id             Number of obs = 611     
 Time variable: year            Min obs per group: 5    
 Number of instruments = 19     Max obs per group: 7    
 Number of groups = 140         Avg obs per group: 5.36 
+------+------------+---------------------+------------+-----------+-----+
|  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
+------+------------+---------------------+------------+-----------+-----+
| L1.n | 0.3912720  |      0.2046422      | 1.9119813  | 0.0558786 |     |
| L2.n | -0.0895707 |      0.0715597      | -1.2516912 | 0.2106824 |     |
|  w   | -0.4985842 |      0.1298172      | -3.8406628 | 0.0001227 | *** |
|  k   | 0.4358161  |      0.0722603      | 6.0311954  | 0.0000000 | *** |
+------+------------+---------------------+------------+-----------+-----+
Hansen test of overid. restrictions: chi(15) = 35.412 Prob > Chi2 = 0.002
Arellano-Bond test for AR(1) in first differences: z = -1.47 Pr > z =0.142
Arellano-Bond test for AR(2) in first differences: z = -0.03 Pr > z =0.977

Dynamic panel-data estimation, two-step system GMM
 Group variable: id             Number of obs = 751     
 Time variable: year            Min obs per group: 5    
 Number of instruments = 27     Max obs per group: 7    
 Number of groups = 140         Avg obs per group: 5.36 
+------+------------+---------------------+------------+-----------+-----+
|  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
+------+------------+---------------------+------------+-----------+-----+
| L1.n | 0.9234415  |      0.2742355      | 3.3673306  | 0.0007590 | *** |
| L2.n | -0.1542449 |      0.1163042      | -1.3262186 | 0.1847673 |     |
|  w   | -0.1867863 |      0.1386538      | -1.3471420 | 0.1779345 |     |
|  k   | 0.1859660  |      0.1350758      | 1.3767531  | 0.1685886 |     |
| _con | 0.8697910  |      0.6687335      | 1.3006543  | 0.1933768 |     |
+------+------------+---------------------+------------+-----------+-----+
Hansen test of overid. restrictions: chi(22) = 55.450 Prob > Chi2 = 0.000
Arellano-Bond test for AR(1) in first differences: z = -1.90 Pr > z =0.058
Arellano-Bond test for AR(2) in first differences: z = -0.30 Pr > z =0.768
4.840353965759277

Their estimates and running times (i.e., total running time of 100 tests) are shown in table below. Scripts of this test are included in the "Benchmark" folder. 

