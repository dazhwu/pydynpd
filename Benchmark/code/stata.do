timer clear
timer on 1
foreach n of numlist 1/100{
clear
insheet using "C:\Users\Tiger\OneDrive\Dynamic Panel\data.csv"
  xtset(id year)
xtabond2 n L(1/2).n w k , gmm(n, lag(2 4)) gmm(w, lag(1 3)) iv(k )  twostep robust 
}
timer off 1

qui timer list
di in r "First time: " r(t1) 
