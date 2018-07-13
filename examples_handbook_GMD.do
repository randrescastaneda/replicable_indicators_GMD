/*==================================================
project:       Examples of the handbook using GMD
Author:        Andres Castaneda 
Dependencies:  The World Bank
----------------------------------------------------
Creation Date:     6 Jul 2018 - 12:03
Modification Date:   
Do-file version:    01
References:          
Output:        http://eca/povdata/examples_handbook_GMD.DO
==================================================*/

/*==================================================
Program set up
==================================================*/
version 15.1
drop _all

* Load data
datalibweb, country(PRY) year(2006 2016) type(GMD) module(ALL)/* 
 */ incppp(welfare) ppp(2011) clear


/*==================================================
deflate welfare variables             
==================================================*/

/* 
There are several welfare variables in the GMD collection. 
In general, the one used for monitoring poverty and inequality is 'welfare'
*/

gen double welfare_ppp = welfare/cpi2011/icp2011/365

table year [w = weight_p], c(mean welfare mean welfare_ppp) format(%11.2f)


/*==================================================
download commands from SSC
==================================================*/

local cmds apoverty ainequal fastgini quantiles prosperity hoi /* 
*/   iop drdecomp skdecomp adecomp

foreach cmd of local cmds {
	capture which `cmd'
	if (_rc != 0) ssc install `cmd'
	else {
		adoupdate `cmd', ssconly
		if ("`r(pkglist)'" == "`cmd'") adoupdate `cmd', update
	}
}

foreach cmd of local cmds {
	capture ssc install `cmd', replace
}

/*==================================================
poverty measures             
==================================================*/

local plines "1.9 3.2 5.5"
local wvar "welfare_ppp"
foreach pl of local plines	{
	gen pl_`=100*`pl'' = `pl'
	forval a=0/2	{
		gen fgt`a'_`=100*`pl'' = 100*((`wvar'<`pl')*(1-(`wvar'/`pl'))^`a')
	}
}

tabstat fgt0* [`w' = weight_p], by(year) nototal
* tabstatmat 
table year [`w' = weight_p], c(mean fgt0_190 mean fgt0_320 mean fgt0_550 )

local w: char _dta[weighttype]
local w = lower("`w'") 

apoverty welfare_ppp [pw = weight_p] if year == 2006, line(1.9)
apoverty welfare_ppp [`w' = weight_p] if year == 2006, varpl(pl_190) 
apoverty welfare_ppp [`w' = weight_p] if year == 2006, varpl(pl_190) all

/*==================================================
Inequality measures              
==================================================*/

ainequal welfare_ppp [w = weight_p] if year == 2006
ainequal welfare_ppp [w = weight_p] if year == 2016, all

/*==================================================
Growth Incidence Curve
==================================================*/
* check command gicurve by Lokshin
bysort year: quantiles welfare_ppp [w = weight_p], /* 
 */          keeptog(hhid) n(100) gen(q)

preserve
rename welfare_ppp* w*

sum year, meanonly
local year1 = r(min)  
local year2 = r(max)  
local period = r(max) - r(min)

sum w [w = weight_p] if year == `year1' , mean
local w1_mean = r(mean)

sum w [w = weight_p] if year == `year2' , mean
local w2_mean = r(mean)

collapse (mean) w [w = weight_p], by(year q)

reshape wide w, i(q) j(year)
gen ann_gic = ((w`year2'/w`year1')^(1/`period')-1)*100
gen ann_gro = ((`w2_mean'/`w1_mean')^(1/`period')-1)*100
sum ann_gic, meanonly
gen avg_gic = r(mean)

label var ann_gic "GIC"
label var ann_gro "Annualized mean growth"
label var avg_gic "Avg. growth of quantiles"

twoway (mspline ann_gic  q, bands(100)) /* 
*/     (line ann_gro avg_gic q), /* 
*/         legend( order(2 3) ) /* 
*/         xtitle("Quantiles") title("Growth Incidence Curve") /* 
*/         subtitle("Paraguay `year1'-`year2'")
restore

/*==================================================
Shared Prosperity
==================================================*/

sum year, meanonly
local year1 = r(min)  
local year2 = r(max)  
local period = r(max) - r(min)

local y = 0
foreach year in `year1' `year2' {
	local ++y
	foreach case in 1 2 3 {
		if      (`case' == 1) local iff ""
		else if (`case' == 2) local iff "& q <= 40"
		else                  local iff "& q  > 40"
		
		sum welfare_ppp [w = weight_p] if year == `year' `iff', meanonly
		local m`y'`case' = r(mean)
		
	}
}

local tot_gr = 100*( (`m21'/`m11')^(1/`period') -1)
local b40_gr = 100*( (`m22'/`m12')^(1/`period') -1)
local t60_gr = 100*( (`m23'/`m13')^(1/`period') -1)
local premium = `b40_gr'-`tot_gr'


display _newline(2) /* 
*/    in g " Annualized overall growth"  _col(32) "= " in w  %5.2f `tot_gr' /* 
*/ _n in g " Annualized growth of B40"   _col(32) "= " in w  %5.2f `b40_gr' /* 
*/ _n in g " Annualized growth of T60"   _col(32) "= " in w  %5.2f `t60_gr' /* 
*/ _n in g " Premium"                    _col(32) "= " in w  %5.2f `premium' _n

prosperity welfare_ppp [w = weight_p], period(year)

/*==================================================
Datt-Ravallion decomposition
==================================================*/

drdecomp welfare_ppp [w = weight_p], by(year) varpl(pl_320)

tempname M
local plines "190 320 550"
qui foreach pl of local plines	{
	drdecomp welfare_ppp [w = weight_p], by(year) /* 
	*/      varpl(pl_`pl') ind(fgt0 fgt1 fgt2)
	mat b = r(b)
	mat `M' = nullmat(`M') \ (J(rowsof(b),1,`pl'), b) // Append results
}

mat colnames `M' = poverty_line fgt effect rate 
svmat double `M', n(col) //converts matrix to variables

* labels for different measures

label define fgt 0 "FGT0" 1 "FGT1" 2 "FGT2", modify
label values fgt fgt

* poverty line labels
label define poverty_line /* 
 */  190 "US $1.9 a day" /* 
 */  320 "US $3.2 a day" /* 
 */  550 "US $5.5 a day" , modify 
label values poverty_line poverty_line

* indicator labels

label define effect 1 "Growth" 2 "Redistribution" 3 "Total change" , modify
label values effect effect

tabdisp poverty_line effect if rate <. , cell(rate) by(fgt) concise // show results
drop poverty_line fgt effect rate 

/*==================================================
Human Opportunity Index
==================================================*/

local vars water banio cloacas electricity landphone cellphone computer /* 
*/        improved_water water_source water_original watertype_quest   /* 
*/        pipedwater_acc improved_sanitation sanitation_source         /* 
*/        sanitation_original toilet_acc

sort hhid
foreach var of local vars {
	cap confirm var `var', exact
	if (_rc) continue
	tempvar m`var'
	by hhid: egen  `m`var'' = max(`var')	
	cap assert `var' == `m`var''
	if (_rc) replace `var' = `m`var''
} 

hoi water            /* 
*/  urban male      /* 
*/  [w=weight_p] if inrange(age, 10, 15) , by(year) decomp2



iop water            /* 
*/  urban male      /* 
*/  [w=weight_p] if inrange(age, 10, 15) & year == 2006, type(d)



exit
/* End of do-file */

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

Notes:
1.
2.
3.


Version Control:


