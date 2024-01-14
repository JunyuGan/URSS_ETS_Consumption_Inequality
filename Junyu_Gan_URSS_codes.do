*******Preamble:* Purpose: This code contains a series of STATA commands designed for empirical analysis of the research topic using DID.* Files Used: [DATA(with co2).dta]* Outputs: The code will be divided into 9 steps: 1. Replace missing value; 2. Winsorization; 3. Construct DID interaction term; 4. Summary statistics for all variables; 5. Baseline regression; 6. Parallel trend test; 7. Robustness test; 8. Heterogeity test; 9. Further discussion* Author: [Junyu Gan]* Date: [2023.8.20]*******

// Set paths
set matsize 11000

// Install packages
capture ssc install winsor
capture ssc install reghdfe
capture ssc install st0085_2 
capture ssc install psmatch2

use "/Users/casillas/Desktop/Distributional effect of energy crisis/DATA (with co2).dta",clear


********************************************************************************

*Step1: Replace missing value

xtset id year
bys id :replace lngdppc=2*l.lngdppc-l2.lngdppc if lngdppc==. //2×previous observation−observation before previous. The expression 2*l.lngdppc-l2.lngdppc effectively calculates the first difference between the two most recent non-missing observations and then adds this difference to the most recent non-missing observation. This can be thought of as a linear extrapolation based on the most recent trend in the data. This is based on the assumption that the most recent trend in the data (as captured by the difference between the two most recent non-missing observations) will continue
bys id :replace lngdppc=2*f.lngdppc-f2.lngdppc if lngdppc==.

bys id :replace uepclf=2*l.uepclf-l2.uepclf if uepclf==.
bys id :replace uepclf=2*f.uepclf-f2.uepclf if uepclf==.

bys id :replace tcpi=2*l.tcpi-l2.tcpi if tcpi==.
bys id :replace tcpi=2*f.tcpi-f2.tcpi if tcpi==.

bys id :replace gc=2*l.gc-l2.gc if gc==.
bys id :replace gc=2*f.gc-f2.gc if gc==.

bys id :replace urbanpop=2*l.urbanpop-l2.urbanpop if urbanpop==.
bys id :replace urbanpop=2*f.urbanpop-f2.urbanpop if urbanpop==.

bys id :replace dependratio=2*l.dependratio-l2.dependratio if dependratio==.
bys id :replace dependratio=2*f.dependratio-f2.dependratio if dependratio==.

********************************************************************************


*Step2: Winsorization

local vv "lngdppc uepclf tcpi gc urbanpop dependratio"
     foreach v of varlist `vv'{
       local a: var lab `v'
	   winsor `v', p(0.05) gen(`v'_x) //winsorization  performed at the 5th percentile from both the lower and upper ends
	   drop `v'
	   rename `v'_x `v'
	   label var `v' "`a'"
     }

********************************************************************************


*Step3：Construct DID interaction term

g dt=(year>=2005)
g dd=post*treatment //post is a dummy variable that takes value of 1 for years >= 2005, 0 otherwise; treatment is a dummy variable that takes value of 1 if the country is included in EU ETS in phase1, 0 otherwise

//control variables
global xlist "lngdppc uepclf tcpi gc urbanpop dependratio"  

******************************************************************************** 


*Step4：Summary statistics for all variables
winsor gini, p(0.05) gen(ginixx)
sum gini dd lngdppc uepclf tcpi gc urbanpop dependratio

****************************************************begin***********************


*Step5：Baseline regression

drop if year>=2008|year<1998  //time frame only up to phase1, phase2 not included
reghdfe ginixx dd, a(id year) cluster(id) 
est store m1
reghdfe ginixx dd $xlist, a(id year) cluster(id)
est store m2
// export results
local model "m1 m2"
esttab `model', mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)   ///
                  scalar(r2 N)     ///
				  compress nogaps	

local model "m1 m2"
esttab `model' using baselineregression.rtf, mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)      ///
                  scalar(r2 N)     ///
				  compress nogaps

********************************************************************************


*Step6：Parallel trend test

//graphical exploration (Figure 1): trends of gini index for treatment and control groups
egen meanp=mean(ginixx), by(year treat)
sort year
twoway (connect meanp year if treatment == 1,lpattern(solid)  lwidth(medthin))   ///
(connect meanp year if treatment == 0,lpattern(longdash) lwidth(medthin)),   ///
 ytitle("ginixx", m(medsmall))   ///
 xline(2005,lpattern(longdash) lwidth(medium)) ///
 xtitle("year", m(medsmall)) legend(label(1 "treatment") label( 2 "control"))  ///
 xscale(range (1998 2007)) ylabel(0.2(0.1)0.4,labs(small)) xlabel(1998 (1) 2007 ,labs(small)) scheme(s1mono)

//event study analysis
gen birth=0 
replace birth = 2005 if treatment==1
gen Dyear=year-birth

//using forvalues loop to generate the "Before" dummy variables and labels
forvalues i=-7(1)-1 {
    local abs_i = abs(`i')
    gen Before`abs_i'=(Dyear==`i' & treatment==1)
    lab var Before`abs_i' "`i'"
}


//generating the "Current" dummy variable and label
gen Current=(Dyear==0 & treatment==1)
lab var Current "0"

//using forvalues loop to generate the "After" dummy variables and labels
forvalues i=1/2 {
    gen After`i'=(Dyear==`i' & treatment==1)
    lab var After`i' "`i'"
}

xtset id year
reghdfe ginixx Before6 Before5 Before4 Before3 Before2 Before1 Current After* $xlist, a(id year) cluster(id)
est store Dynamic

local model "Dynamic"
esttab `model' using paralleltrend.rtf, mtitle(`model') replace ///
      b(%7.3f) se(%7.3f) ///
      star(* 0.1 ** 0.05 *** 0.01) ///
      scalar(r2 N) ///
      compress nogaps


// Figure 2
coefplot Dynamic, ///
keep(Before* Current After*)   ///
vertical  ///
 yline(0, lwidth(vthin) lpattern(dash) lcolor(teal)) ///
 xline(0, lwidth(vthin) lpattern(dash) lcolor(teal)) ///
   ylabel(-0.01(0.002)0.002) ///
   xline(0, lwidth(vthin) lpattern(dash) lcolor(teal)) ///
   ytitle("consumption gini", size(small)) ///
   xtitle("year", size(small)) ///
   addplot(line @b @at , lcolor(gs0) lwidth(thin))   ///
   ciopts(lpattern(dash) recast(rcap) msize(medium))  ///
   msymbol(circle_hollow) ///
   scheme(s1mono) level(90)

********************************************************************************


*Step7；Robustness test

**(1)Expected effect test (exclude pre treatment impact of the policy)
xtset id year
set matsize 11000
g d2004=(year==2004)
reghdfe ginixx dd i.treatment#i.d2004 $xlist, a(id year) cluster(id)  //coefficient of `i.treatment#i.d2004` check for any differential change in ginixx for treated entities in the year 2004. Adding it to the regression to control for the pre treatment effect of the policy.
est store m1
local model "m1"
esttab `model', mtitle(`model') replace   ///
	  b(%7.3f) se(%6.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)    ///
                  scalar(r2 N N_g)     ///
				  compress nogaps
local model "m1"
esttab `model' using expectedeffect.rtf, mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)      ///
                  scalar(r2 N)     ///
				  compress nogaps 
********************************************************************************

**(2)PSM-did
global ylist "lngdppc uepclf tcpi gc urbanpop dependratio"   
set seed 0001	//define seeds
gen tmp = runiform() //generate random number
sort tmp //sorts the dataset based on the random number,
psmatch2 treatment $ylist, out(ginixx) logit ate kernel k(biweight) common caliper(.05) ties //matches treated observations with untreated ones based on the propensity scores
pstest $ylist, both graph  //tests if the covariates in $ylist are balanced between the treated and untreated groups after matching. 
gen common=_support
psgraph,scheme(s1mono)
xtset id year
reghdfe ginixx dd  $xlist if common==1, a(id year) cluster(id)
est store m
local model "m"
esttab `model', mtitle(`model') replace   ///
	  b(%7.3f) se(%6.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)   ///
                  scalar(r2 N N_g)     ///
				  compress nogaps
local model "m"
esttab `model' using PSM-dd.rtf, mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)      ///
                  scalar(r2 N)     ///
				  compress nogaps 
*******************************************************************************

**(3）Alternative index
winsor palmaratio, p(0.05) gen(palmaratiox) //use Palmaratio instead of gini index
reghdfe palmaratiox dd $xlist, a(id year) cluster(id)
est store m

local model "m"
esttab `model', mtitle(`model') replace   ///
	  b(%7.3f) se(%6.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)   ///
                  scalar(r2 N N_g)     ///
				  compress nogaps
local model "m"
esttab `model' using Alternativeindex.rtf, mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)    ///
                  scalar(r2 N)     ///
				  compress nogaps 
********************************************************************************


*Step8: Heterogeneity test

sum lngdppc if year == 2004
local mean_lngdppc = r(mean)
gen quart_lngdppc_2004 = (lngdppc > `mean_lngdppc' & year==2004) //construct only for the year you take the mean in
bys id: egen quart_lngdppc=max(quart_lngdppc_2004) //take maximum by country -- this ensures a country is either 0 or 1 for all years
drop quart_lngdppc_2004
gen lngdppcdd=dd*quart_lngdppc
reghdfe ginixx lngdppcdd dd uepclf tcpi gc urbanpop dependratio, a(id year) cluster(id)
est store m1

sum ecpi if year == 2004
local mean_ecpi = r(mean)
gen quart_ecpi_2004 = (ecpi > `mean_ecpi' & year==2004)
bys id: egen quart_ecpi=max(quart_ecpi_2004)
drop quart_ecpi_2004
gen ecpidd=dd*quart_ecpi
reghdfe ginixx ecpidd dd lngdppc uepclf tcpi gc urbanpop dependratio, a(id year) cluster(id)
est store m2

sum lnco2 if year == 2004
local mean_lnco2 = r(mean)
gen quart_lnco2_2004 = (lnco2 > `mean_lnco2' & year==2004)
bys id: egen quart_lnco2=max(quart_lnco2_2004)
drop quart_lnco2_2004
gen lnco2dd=dd*quart_lnco2
reghdfe ginixx lnco2dd dd lngdppc uepclf tcpi gc urbanpop dependratio, a(id year) cluster(id)
est store m3

sum tcpi if year == 2004
local mean_tcpi = r(mean)
gen quart_tcpi_2004 = (tcpi > `mean_tcpi' & year==2004)
bys id: egen quart_tcpi=max(quart_tcpi_2004)
drop quart_tcpi_2004
gen tcpidd=dd*quart_tcpi
reghdfe ginixx tcpidd dd lngdppc uepclf gc urbanpop, dependratio, a(id year) cluster(id)
est store m4

local model "m1 m2 m3 m4"
esttab `model', mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)    ///
                  scalar(r2 N)     ///
				  compress nogaps
				  
local model "m1 m2 m3 m4"
esttab `model' using 异质性1.rtf, mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)        ///
                  scalar(r2 N)     ///
				  compress nogaps 
 
********************************************************************************


*Step9: Further discussion (Include phase2 countries and years)

use "C:\Users\28124\Desktop\further discussion.dta", clear //interaction term dd2 already defined in dataset. dd2 takes the value of 1 only when: 1. the country is included in the EU ETS & 2. the year of the observation is post treatment
xtset id year
bys id :replace lngdppc=2*l.lngdppc-l2.lngdppc if lngdppc==.
bys id :replace lngdppc=2*f.lngdppc-f2.lngdppc if lngdppc==.

bys id :replace uepclf=2*l.uepclf-l2.uepclf if uepclf==.
bys id :replace uepclf=2*f.uepclf-f2.uepclf if uepclf==.

bys id :replace tcpi=2*l.tcpi-l2.tcpi if tcpi==.
bys id :replace tcpi=2*f.tcpi-f2.tcpi if tcpi==.

bys id :replace gc=2*l.gc-l2.gc if gc==.
bys id :replace gc=2*f.gc-f2.gc if gc==.

bys id :replace urbanpop=2*l.urbanpop-l2.urbanpop if urbanpop==.
bys id :replace urbanpop=2*f.urbanpop-f2.urbanpop if urbanpop==.

bys id :replace dependratio=2*l.dependratio-l2.dependratio if dependratio==.
bys id :replace dependratio=2*f.dependratio-f2.dependratio if dependratio==.

local vv "lngdppc uepclf tcpi gc urbanpop dependratio"
     foreach v of varlist `vv'{
       local a: var lab `v'
	   winsor `v', p(0.05) gen(`v'_x)
	   drop `v'
	   rename `v'_x `v'
	   label var `v' "`a'"
     }

winsor gini, p(0.05) gen(ginixx)

sum gini dd lngdppc uepclf tcpi gc urbanpop dependratio
reghdfe ginixx dd2 if year<2013&year>=1998, a(id year) cluster(id) 
est store m1
reghdfe ginixx dd2 lngdppc uepclf tcpi gc urbanpop dependrati if year<2013&year>=1998, a(id year) cluster(id)
est store m2

//export result
local model "m1 m2"
esttab `model', mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)   ///
                  scalar(r2 N)     ///
				  compress nogaps	

local model "m1 m2"
esttab `model' using further discussion.rtf, mtitle(`model') replace   ///
	  b(%7.3f) se(%7.3f)   ///
                  star(* 0.1 ** 0.05 *** 0.01)      ///
                  scalar(r2 N)     ///
				  compress nogaps





