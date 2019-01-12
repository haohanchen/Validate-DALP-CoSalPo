******************************************STATA CODE FOR 2010 APSA PAPER ON MEASURES OF PROGRAMMATIC PARTY COMPETITION****************************************
*************************************************************KENT FREEZE AND HERBERT KITSCHELT****************************************************************

*QUESTIONS OR COMMENTS ON THIS DO-FILE SHOULD BE DIRECTED TO KENT FREEZE AT FREEZEKENT@GMAIL.COM OR KEF13@DUKE.EDU

*RUN USING STATA 11.1.  SOME CODE MAY NOT WORK FOR VERSIONS OF STATA PREVIOUS TO 11.0

*This do-file uses the user-written package "ROWRANKS", ensure that it is installed.

ssc install rowranks

*Prior to running this do-file, change the working directory to your local drive!!!!!!


cd "C:\Dropbox\Data\democraticaccountability"

use "cumulative_36partysize.dta", clear


********************************************Number of Issue Categories and number of parties Per Country*************************************************

aorder

foreach x of varlist  d7_ago-d21_vene {
rename `x' dz_`x'
}
aorder 
order  country cname ccode ccodewb date time authenticateduserid journalist PARTY party
drop d22-dy

collapse (mean) d1-dz_d21_vene, by(country cname ccode ccodewb PARTY party)

*Generation of the Number of issues in each country

	gen issuenum=0
	foreach x of varlist d1- dz_d21_vene {
	replace issuenum=issuenum+1 if `x'!=.
	}
	
	drop d1- dz_d21_vene
	
	collapse (max) issuenum, by(cname)


merge 1:m cname using "cumulative_36partysize.dta"
drop _merge
drop if party==77

egen partynum=max(party), by(ccode)

***********************************************************************************************************************************************
*****************************************************PARTY-LEVEL CALCULATIONS******************************************************************
***********************************************************************************************************************************************


						*****************Individual-Level Calculations to Party-Level Variables*********************
						
******************************************SECTION A: SETUPS************************************************************************************

*Getting D-section variables in order for program calculations

aorder

foreach x of varlist  d7_ago-d21_vene {
rename `x' dz_`x'
}

drop d22-d46
renpfix dz_

aorder 

order  country cname ccode ccodewb date time authenticateduserid journalist PARTY party

*Need to shrink the mess of d variables into something more manageable, otherwise the calculations tend to run slowly and potentially out of space

gen sum=0

foreach w in 6 7 8 9 10 11 12 13 14 15 16 17 18 19 {
gen d`w'=.
}

foreach x of varlist d7_ago-d21_vene {
egen mean`x'= mean(`x'), by(ccode)
replace mean`x'=1 if mean`x'!=.
foreach w in 6 7 8 9 10 11 12 13 14 15 16 17 18 19 {
replace sum=1 if mean`x'==1 & d`w'==.
replace d`w'=`x' if mean`x'==1 & d`w'==.
replace mean`x'=. if sum==1
replace sum=0
}
}

order d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17 d18 d19

keep country cname ccode ccodewb party PARTY authenticateduserid d1-d19 issuenum partynum p31

aorder


**************************************************Correction for personal-anchor points*******************************************

*Remove the mean of each expert's judgments for the parties she judged.
*Transform expert judgments into standardized z-scores, factoring in the extremism of certain experts.
	
	foreach x of varlist d1-d19 {
	gen original`x'=`x'
	drop `x' 
	gen `x'=original`x' /*This was just to quickly get rid of annoying labels for d questions*/
	replace `x'=. if `x'==88
	replace `x'=. if `x'==99
	egen expmean`x'=mean(`x'), by(ccode authenticateduserid)
	egen expsd`x'=sd(`x'), by(ccode authenticateduserid)
	replace `x'=(`x'-expmean`x')/expsd`x'
	drop original`x'
	drop expmean`x'
	}
	
	order d1-d19
	egen avgscore=rowmean(d1-d19)
	drop if avgscore==. /*Getting rid of responses for which all values are missing*/
	drop avgscore
	aorder

	
****************************************************SECTION B: Running of Discriminant Analysis***************************************************

*STEP 1: Imputation of Missing Data
*For some countries, the number of valid responses is insufficient to run the analysis.  So for the discriminant analysis (and only the
*discriminant analysis, we used imputed data.

*This is a "quick and dirty" way to impute data - however Stata's multiple impute command does not work well with "hierarchical" data
*such as is being used here.  In the future, this would be a good thing to correct.  But since the imputation isn't really used to
*test any theories, but just to build discriminant functions, it is less of a concern.

foreach y of numlist 1 2 to 19 { /*Copying original variables for imputation*/
gen imputed`y'=d`y'
egen meand`y'=mean(d`y'), by(ccode party) /*Generating Means for each party*/
egen standdevd`y'=sd(d`y'), by(ccode party) /*Generating SD for each party*/
replace imputed`y'=rnormal(meand`y', standdevd`y') if imputed`y'==. & `y'<(issuenum+1) /*Imputing using Normal Dist. and Mean and SD for each party*/
}


/*STEP 2: RUNNING DISCRIMINANT ANALYSIS

*NOTE: HONDURAS (CCODE 340) DID NOT HAVE ENOUGH RESPONSES, EVEN AFTER IMPUTATION, LARGELY DUE TO THE LARGE NUMBER OF ISSUES THERE (19).

gen ccodes=ccode
replace ccodes=. if ccode==158

aorder

gen cancor=.
gen eigen=.
forvalues x = 1/19 {
gen dscldd`x'=.
}

levelsof ccodes, local(ccodes)
foreach x of local ccodes {
display `x'
quietly sum issuenum if ccode==`x'
local isnum=r(mean)
candisc imputed1-imputed`isnum' if ccode==`x', group(party)
matrix discrim=e(candisc_stat)
svmat discrim, names( col )
egen maxcancor=max(Canonical_cor)
egen maxeigen=max(Eigenvalue)

replace cancor=maxcancor if ccode==`x' /*For each country return the strongest factor's canonical correlation and eigenvalue*/
replace eigen=maxeigen if ccode==`x'

matrix canstr=e(canstruct) /*Bringing in the Issue discriminant loadings*/
svmat canstr, names( col )

local funcnum=e(f) /*need a macro for the number of functions*/

forvalues y = 1/`funcnum' { /*marking functions for which the pvalue is less than 0.05 as missing*/
egen pvaluefunc`y'=max(pvalue) if _n==`y'
egen pvalue2func`y'=max(pvaluefunc`y')
replace function`y'=. if pvalue2func`y'>0.05
replace function`y'=(-1)*function`y' if function`y'<0
}
egen maxdiscloadings=rowmax(function1-function`funcnum')

forvalues z = 1/`isnum' { /*putting the maximum discriminant loading for each issue into the correct variable*/
egen discloading`z'=max(maxdiscloadings) if _n==`z'
egen disc2loading`z'=max(discloading`z')
replace dscldd`z'=disc2loading`z' if ccode==`x'
replace dscldd`z'=0 if dscldd`z'==. & ccode==`x'
}

drop  Canonical_corr-disc2loading`isnum'

}

drop  imputed1- imputed19  meand1- meand19 standdevd1- standdevd19
*/

**************************************SECTION C: CREATION OF COSALPO SCORES************************************************************

*1. Generating Salience Measure
	*Variable=% of experts who rated a party on an issue.

	* Generate percent nonmissing for each of the policy positions
	
	gen constant=1 if party
	egen constantct=count(constant), by(ccode party)
	foreach x of varlist d1-d19 {
	egen ct`x'=count(`x'), by(ccode party)
	gen nr`x'=ct`x'/constantct
	}
	
	aorder
	
	drop ctd1-ctd19
	drop constant
	drop constantct


*2. Generating Cohesion Measure
	*Variable = sd of expert judgments on a particular issue.

	* Generate standard deviations for each of the policy positions
	foreach x of varlist d1-d19 {
	egen sd`x'=sd(`x'), by(ccode party)
	}
	
	aorder

***********************************************Calculations at the Party-Level*******************************************

collapse (mean) issuenum sdd1-sdd19 nrd1-nrd19 p31 d1-d19 partynum, by(country cname ccode ccodewb PARTY party)
aorder

*NOTE: THE PARTY LEVEL POLARIZATION MEASURE IS CALCULATED BY CONVERTING THE DATA TO WIDE FORM, WHICH CAUSES THE DATA TO RUN OUT OF 
*CALCULATION SPACE WITH TOO MANY VARIABLES.  THEREFORE, ONLY ESSENTIAL VARIABLES ARE KEPT, WITH OTHER VARIABLES SAVED SEPARATELY AND
*LATER MERGED BACK IN TO THE DATA SET.

save "tempcalc.dta", replace

use "tempcalc.dta", clear

keep ccode party p31 d1-d19

*3. Generating Programmatic Polarization Variable
	
	order ccode party p31
	foreach x of varlist ccode-d19 {
	rename `x' `x'_
	}

	reshape wide p31_ d1_-d19_, i(ccode_) j(party)

	aorder
	
	foreach x of varlist p31_1-p31_17 {
	foreach y of varlist p31_1-p31_17 {
	gen `x'`y'sum=`x' + `y'
	}
	replace `x'`x'sum=.
	egen `x'tot = rowtotal(`x'p31_1sum-`x'p31_17sum)
	}
	
	foreach w of numlist 1 2 to 19 {
	foreach x of numlist 1 2 to 17 {
	foreach y of numlist 1 2 to 17 {
	gen d`w'_`x'd`w'_`y'diff= d`w'_`x'-d`w'_`y' if d`w'_`x'!=d`w'_`y'
	replace d`w'_`x'd`w'_`y'diff=(-1)*d`w'_`x'd`w'_`y'diff if d`w'_`x'd`w'_`y'diff<0
	replace d`w'_`x'd`w'_`y'diff=d`w'_`x'd`w'_`y'diff*p31_`x'p31_`y'sum
	}
	egen partypolard`w'_`x'= rowtotal(d`w'_`x'd`w'_1diff-d`w'_`x'd`w'_17diff)
	replace partypolard`w'_`x'=partypolard`w'_`x'/p31_`x'tot
	order partypolard`w'_`x'
	}
	drop  d`w'_1d`w'_1diff- d`w'_17d`w'_17diff
	}
		
	drop p31_1p31_1sum- p31_17tot
	
	order _all, sequential
	reshape long 
	drop if ccode==.
	
	rensfix _
	
	foreach x of numlist 1 2 to 19 {
	gen parpold`x'=.
	foreach y of numlist 1 2 to 17 {
	replace parpold`x'= partypolard`x'_`y' if party==`y'
	}
	}
	
	drop partypolard1_1- partypolard19_17
	keep ccode party parpold1- parpold19
	drop if parpold1==.
	
	merge 1:1 ccode party using "tempcalc.dta"
	aorder
	drop _merge

**********************************************************Normalizations***************************************************************

*1. Normalization for Cohesion Measure.

aorder

egen sdmax=rowmean(sdd1-sdd19)
egen sdsupermax=max(sdmax)
egen sdlowmax=min(sdmax)

foreach x of varlist sdd1-sdd19 {
replace `x'=0 if `x'<=sdlowmax
replace `x'=`x'/sdsupermax        		
replace `x'=1 if `x'>1
replace `x'=1-`x'        			/*Reverse order so greater values reflect greater Cohesion*/
}

*2. Normalization for Polarization Measure

egen parpolmax=rowmean(parpold1-parpold19)
egen parpolsupermax=max(parpolmax)
egen parpollowmax=min(parpolmax)

foreach x of varlist parpold1-parpold19 {
replace `x'=0 if `x'<=parpollowmax
replace `x'=`x'/parpolsupermax
replace `x'=1 if `x'>1
}

*3. Normalization for Salience Measure

foreach x of varlist nrd1-nrd19 {
replace `x'=(`x'-0.4)/0.6
replace `x'=0 if `x'<0
}

*Making it so that nr parpol and sd measures are blank when they should be
foreach x of numlist 1 2 to 19 {
replace nrd`x'=. if d`x'==.
replace sdd`x'=. if d`x'==.
replace parpold`x'=. if d`x'==.
}

*****************************************************GENERATION OF PPS SCORES************************************************************

*1. CoSalPo Score
foreach x of numlist 1 2 to 19 {
gen coposal_d`x'=nrd`x'*sdd`x'*parpold`x'
}

*2. CoSal Score
foreach x of numlist 1 2 to 19 {
gen cosal_d`x'=nrd`x'*sdd`x'
}

aorder

*******************************************************SECTION D: AGGREGATIONS TO THE PARTY-LEVEL***************************************************

rowranks (coposal_d1-coposal_d19), gen(coposal_rankd1-coposal_rankd19) field
rowranks (cosal_d1-cosal_d19), gen(cosal_rankd1-cosal_rankd19) field

*1. Average PPS over all issues in the country
	*a) Coposal Score
	egen copolsal_mean=rowmean(coposal_d1-coposal_d19)
	*b) Cosal Score
	egen cosal_mean=rowmean(cosal_d1-cosal_d19)
	
*2. 8-Issues: 
	*i) D1 or D3, but not both.
	*ii) D4
	*iii) D5
	*iv) Up to 5 more issues with the highest PPS scores.

	*a) Coposal Score
	foreach x of numlist 1 2 to 8 {
	gen coposal_8_`x'=.
	}
	
	replace coposal_8_1=coposal_d1 if coposal_rankd1>=coposal_rankd3
	replace coposal_8_1=coposal_d3 if coposal_rankd1<coposal_rankd3
	replace coposal_8_2=coposal_d4
	replace coposal_8_3=coposal_d5
	rowranks (coposal_d6-coposal_d19), gen(coposal_rank8d6-coposal_rank8d19) field
	
	foreach x of numlist 6 7 to 19 {
	replace coposal_8_4=coposal_d`x' if coposal_rank8d`x'==1
	replace coposal_8_5=coposal_d`x' if coposal_rank8d`x'==2
	replace coposal_8_6=coposal_d`x' if coposal_rank8d`x'==3
	replace coposal_8_7=coposal_d`x' if coposal_rank8d`x'==4
	replace coposal_8_8=coposal_d`x' if coposal_rank8d`x'==5
	}
	
	egen coposal_8=rowmean(coposal_8_1-coposal_8_8)
	
	*b) Cosal Score
	
	foreach x of numlist 1 2 to 8 {
	gen cosal_8_`x'=.
	}
	
	replace cosal_8_1=cosal_d1 if cosal_rankd1>=cosal_rankd3
	replace cosal_8_1=cosal_d3 if cosal_rankd1<cosal_rankd3
	replace cosal_8_2=cosal_d4
	replace cosal_8_3=cosal_d5
	rowranks (cosal_d6-cosal_d19), gen(cosal_rank8d6-cosal_rank8d19) field
	
	foreach x of numlist 6 7 to 19 {
	replace cosal_8_4=cosal_d`x' if cosal_rank8d`x'==1
	replace cosal_8_5=cosal_d`x' if cosal_rank8d`x'==2
	replace cosal_8_6=cosal_d`x' if cosal_rank8d`x'==3
	replace cosal_8_7=cosal_d`x' if cosal_rank8d`x'==4
	replace cosal_8_8=cosal_d`x' if cosal_rank8d`x'==5
	}
	
	egen cosal_8=rowmean(cosal_8_1-cosal_8_8)


*3. 4-Issues:
	*i) D1 and D3 OR D2 and D3, or one of them.
	*ii) D4
	*iii) D5
	*iv) Most-salient country-specific issue if PPS> 2 highest of D1-D3.
	
	*a) Coposal Score
	
	foreach x of numlist 1 2 to 4 {
	gen coposal_4_`x'=.
	}
	
	rowranks (coposal_d1-coposal_d3), gen(coposal_rank4d1-coposal_rank4d3) field
	
	foreach x of numlist 1 2 3 {
	replace coposal_4_1=coposal_d`x' if coposal_rank4d`x'==1
	}
	replace coposal_4_1=coposal_d3 if coposal_4_1==. /* there are a few cases with no rank of 1 among the first three issues - they all tie at `2' */
	replace coposal_4_2=coposal_d4
	replace coposal_4_3=coposal_d5
	
	egen tempblah=rowmin(coposal_rank8d6-coposal_rank8d19)
	
	egen coposal_rank4top=rowmax(coposal_d6-coposal_d19)
	egen tempmax=rowmax(coposal_d1-coposal_d3 coposal_d6-coposal_d19)
	
	foreach x of numlist 6 7 to 19 {
	replace coposal_4_4=coposal_d`x' if coposal_rank8d`x'==tempblah
	}

	foreach x of numlist 1 2 to 3 {
	replace coposal_4_4=coposal_d`x' if coposal_d`x'>coposal_4_4 & coposal_rank4d`x'!=3 & coposal_rank4d`x'!=1 | coposal_4_4==. & coposal_rank4d`x'!=1 & coposal_rank4d`x'!=3
	}
	
	egen coposal_4=rowmean(coposal_4_1-coposal_4_4)

	*b) Cosal Measure:
	
	foreach x of numlist 1 2 to 4 {
	gen cosal_4_`x'=.
	}
	
	rowranks (cosal_d1-cosal_d3), gen(cosal_rank4d1-cosal_rank4d3) field
	
	foreach x of numlist 1 2 3 {
	replace cosal_4_1=cosal_d`x' if cosal_rank4d`x'==1
	}
	replace cosal_4_2=cosal_d4
	replace cosal_4_3=cosal_d5
	
	gen cosal_rank4top=.
	foreach x of numlist 6 7 to 19 {
	replace cosal_rank4top=cosal_d`x' if cosal_rank8d`x'==1
	}
	
	foreach x of numlist 1 2 3 {
	replace cosal_4_4=cosal_d`x' if cosal_rank4d`x'==2
	}
	
	replace cosal_4_4=cosal_rank4top if cosal_rank4top>cosal_4_4
	
	replace cosal_4_4=cosal_rank4top if cosal_rank4d2==3
	
	egen cosal_4=rowmean(cosal_4_1-cosal_4_4)

*4. 3-Issues:
	*i) One of D1-D3
	*ii) D4
	*iii) D5
 
	*a) Coposal Measure:
	
	foreach x of numlist 1 2 3 {
	gen coposal_3same_`x'=.
	replace coposal_3same_1=coposal_d`x' if coposal_rank4d`x'==1
	}
	replace coposal_3same_2=coposal_d4
	replace coposal_3same_3=coposal_d5
	
	egen coposal_3same=rowmean(coposal_3same_1-coposal_3same_3)
	
	*b) Cosal Measure:
	
	foreach x of numlist 1 2 3 {
	gen cosal_3same_`x'=.
	replace cosal_3same_1=cosal_d`x' if cosal_rank4d`x'==1
	}
	replace cosal_3same_2=cosal_d4
	replace cosal_3same_3=cosal_d5
	
	egen cosal_3same=rowmean(cosal_3same_1-cosal_3same_3)


*5. 3-Issues, Modified:
	*i) all of D1-D3.

	*a) Coposal Measure:
	egen coposal_3econ=rowmean(coposal_d1-coposal_d3)

	*b) Cosal Measure:
	egen cosal_3econ=rowmean(cosal_d1-cosal_d3)

order country cname ccode ccodewb party PARTY 

drop  parpollowmax- parpolsupermax  sdlowmax- cosal_rankd19  coposal_8_1- coposal_8_8  cosal_8_1- cosal_8_8  coposal_4_1- coposal_rank4top  cosal_4_1- cosal_rank4top  coposal_3same_1- coposal_3same_3  cosal_3same_1- cosal_3same_3

egen coposalmax=rowmax(coposal_d1-coposal_d17)

egen cosalmax=rowmax(cosal_d1-cosal_d17)

aorder

order country cname ccode ccodewb party PARTY

drop  coposal_rank8d6- coposal_rank8d19  cosal_rank8d6- cosal_rank8d19 tempmax

save "freezekitscheltpartylevel_mar2013.dta", replace


***********************************************************************************************************************************************
*********************************************************NATIONAL-LEVEL CALCULATIONS***********************************************************
***********************************************************************************************************************************************

use "freezekitscheltpartylevel_mar2013.dta", clear

foreach x of varlist d1-d19 {
replace parpol`x'=`x'
}

*Creation of National-Level Data Set.  Note that parpold1-parpold19 have been transformed into the country-level polarization measure.  

collapse (mean) issuenum nrd1- nrd19  sdd1- sdd19 d1-d19 partynum (sd)  parpold1-parpold19  [aweight=p31], by( country cname ccode ccodewb)


						**************Normalization for Polarization Measure********************

egen parpolmax=rowmean(parpold1-parpold19)
egen parpolsupermax=max(parpolmax)
egen parpollowmax=min(parpolmax)

foreach x of varlist parpold1-parpold19 {
replace `x'=0 if `x'<=parpollowmax
replace `x'=`x'/parpolsupermax
replace `x'=1 if `x'>1
}

*Making it so that parpol measures are blank when they should be
foreach x of numlist 1 2 to 19 {
replace parpold`x'=. if d`x'==.
}

*****************************************************GENERATION OF PPS SCORES************************************************************

*1. CoPoSal Score
foreach x of numlist 1 2 to 19 {
gen coposal_d`x'=nrd`x'*sdd`x'*parpold`x'
}

*2. CoSal Score
foreach x of numlist 1 2 to 19 {
gen cosal_d`x'=nrd`x'*sdd`x'
}


aorder

**********************************************AGGREGATIONS TO CREATE COUNTRY-LEVEL INDICATORS*********************************************
rowranks (coposal_d1-coposal_d19), gen(coposal_rankd1-coposal_rankd19) field
rowranks (cosal_d1-cosal_d19), gen(cosal_rankd1-cosal_rankd19) field

*1. Average PPS over all issues in the country
	*a) Coposal Score
	egen copolsal_mean=rowmean(coposal_d1-coposal_d19)
	*b) Cosal Score
	egen cosal_mean=rowmean(cosal_d1-cosal_d19)
	
*2. 8-Issues: 
	*i) D1 or D3, but not both.
	*ii) D4
	*iii) D5
	*iv) Up to 5 more issues with the highest PPS scores.

	*a) Coposal Score
	foreach x of numlist 1 2 to 8 {
	gen coposal_8_`x'=.
	}
	
	replace coposal_8_1=coposal_d1 if coposal_rankd1>=coposal_rankd3
	replace coposal_8_1=coposal_d3 if coposal_rankd1<coposal_rankd3
	replace coposal_8_2=coposal_d4
	replace coposal_8_3=coposal_d5
	rowranks (coposal_d6-coposal_d19), gen(coposal_rank8d6-coposal_rank8d19) field
	

	foreach x of numlist 6 7 to 19 {
	replace coposal_8_4=coposal_d`x' if coposal_rank8d`x'==1
	replace coposal_8_5=coposal_d`x' if coposal_rank8d`x'==2
	replace coposal_8_6=coposal_d`x' if coposal_rank8d`x'==3
	replace coposal_8_7=coposal_d`x' if coposal_rank8d`x'==4
	replace coposal_8_8=coposal_d`x' if coposal_rank8d`x'==5
	}
	
	egen coposal_8=rowmean(coposal_8_1-coposal_8_8)
	
	*b) Cosal Score
	
	foreach x of numlist 1 2 to 8 {
	gen cosal_8_`x'=.
	}
	
	replace cosal_8_1=cosal_d1 if cosal_rankd1>=cosal_rankd3
	replace cosal_8_1=cosal_d3 if cosal_rankd1<cosal_rankd3
	replace cosal_8_2=cosal_d4
	replace cosal_8_3=cosal_d5
	rowranks (cosal_d6-cosal_d19), gen(cosal_rank8d6-cosal_rank8d19) field
	
	foreach x of numlist 6 7 to 19 {
	replace cosal_8_4=cosal_d`x' if cosal_rank8d`x'==1
	replace cosal_8_5=cosal_d`x' if cosal_rank8d`x'==2
	replace cosal_8_6=cosal_d`x' if cosal_rank8d`x'==3
	replace cosal_8_7=cosal_d`x' if cosal_rank8d`x'==4
	replace cosal_8_8=cosal_d`x' if cosal_rank8d`x'==5
	}
	
	egen cosal_8=rowmean(cosal_8_1-cosal_8_8)

*3. 4-Issues:
	*i) D1 and D3 OR D2 and D3, or one of them.
	*ii) D4
	*iii) D5
	*iv) Most-salient country-specific issue if PPS> 2 highest of D1-D3.
	
	*a) Coposal Score
	
	foreach x of numlist 1 2 to 4 {
	gen coposal_4_`x'=.
	}
	
	rowranks (coposal_d1-coposal_d3), gen(coposal_rank4d1-coposal_rank4d3) field
	
	foreach x of numlist 1 2 3 {
	replace coposal_4_1=coposal_d`x' if coposal_rank4d`x'==1
	}
	replace coposal_4_1=coposal_d3 if coposal_4_1==. /* there are a few cases with no rank of 1 among the first three issues - they all tie at `2' */
	replace coposal_4_2=coposal_d4
	replace coposal_4_3=coposal_d5
	
	egen tempblah=rowmin(coposal_rank8d6-coposal_rank8d19)
	
	egen coposal_rank4top=rowmax(coposal_d6-coposal_d19)
	egen tempmax=rowmax(coposal_d1-coposal_d3 coposal_d6-coposal_d19)
	
	foreach x of numlist 6 7 to 19 {
	replace coposal_4_4=coposal_d`x' if coposal_rank8d`x'==tempblah
	}

	foreach x of numlist 1 2 to 3 {
	replace coposal_4_4=coposal_d`x' if coposal_d`x'>coposal_4_4 & coposal_rank4d`x'!=3 & coposal_rank4d`x'!=1 | coposal_4_4==. & coposal_rank4d`x'!=1 & coposal_rank4d`x'!=3
	}
	
	egen coposal_4=rowmean(coposal_4_1-coposal_4_4)

	
	*b) Cosal Measure:
	
	foreach x of numlist 1 2 to 4 {
	gen cosal_4_`x'=.
	}
	
	rowranks (cosal_d1-cosal_d3), gen(cosal_rank4d1-cosal_rank4d3) field
	
	foreach x of numlist 1 2 3 {
	replace cosal_4_1=cosal_d`x' if cosal_rank4d`x'==1
	}
	replace cosal_4_2=cosal_d4
	replace cosal_4_3=cosal_d5
	
	gen cosal_rank4top=.
	foreach x of numlist 6 7 to 19 {
	replace cosal_rank4top=cosal_d`x' if cosal_rank8d`x'==1
	}
	
	foreach x of numlist 1 2 3 {
	replace cosal_4_4=cosal_d`x' if cosal_rank4d`x'==2
	}
	
	replace cosal_4_4=cosal_rank4top if cosal_rank4top>cosal_4_4
	
	replace cosal_4_4=cosal_rank4top if cosal_rank4d2==3
	
	egen cosal_4=rowmean(cosal_4_1-cosal_4_4)
	
*4. 3-Issues:
	*i) One of D1-D3
	*ii) D4
	*iii) D5
 
	*a) Coposal Measure:
	
	foreach x of numlist 1 2 3 {
	gen coposal_3same_`x'=.
	replace coposal_3same_1=coposal_d`x' if coposal_rank4d`x'==1
	}
	replace coposal_3same_2=coposal_d4
	replace coposal_3same_3=coposal_d5
	
	egen coposal_3same=rowmean(coposal_3same_1-coposal_3same_3)
	
	*b) Cosal Measure:
	
	foreach x of numlist 1 2 3 {
	gen cosal_3same_`x'=.
	replace cosal_3same_1=cosal_d`x' if cosal_rank4d`x'==1
	}
	replace cosal_3same_2=cosal_d4
	replace cosal_3same_3=cosal_d5
	
	egen cosal_3same=rowmean(cosal_3same_1-cosal_3same_3)

*5. 3-Issues, Modified:
	*i) all of D1-D3.

	*a) Coposal Measure:
	egen coposal_3econ=rowmean(coposal_d1-coposal_d3)

	*b) Cosal Measure:
	egen cosal_3econ=rowmean(cosal_d1-cosal_d3)

aorder

egen coposalmax=rowmax(coposal_d1-coposal_d17)

egen cosalmax=rowmax(cosal_d1-cosal_d17)

gen maxdnum=""

forvalues x = 1/19 {
replace maxdnum="d`x'" if coposalmax==coposal_d`x'
}

order country cname ccode ccodewb copolsal_mean cosal_mean coposal_8 cosal_8 coposal_4 cosal_4 coposal_3same cosal_3same coposal_3econ cosal_3econ coposalmax cosalmax

drop  coposal_rankd1-cosal_rankd19 parpollowmax- parpolsupermax tempmax



save "freezekitscheltcountrylevel_mar2013.dta", replace

*Calculation Party system Fractionalization
use "freezekitscheltpartylevel_mar2013.dta", replace

gen p31sq=(p31/100)^2

egen partyfract=sum(p31sq), by(ccode)
replace partyfract=1-partyfract

collapse partyfract, by(ccode)

merge 1:1 ccode using "freezekitscheltcountrylevel_mar2013.dta
drop _merge

save "freezekitscheltcountrylevel_mar2013.dta", replace

keep ccode party coposal_4 coposal_3econ 
rename coposal_4 coposal_4_stand
rename coposal_3econ coposal_3econ_stand



/*Generating Country-Level Scattergrams, etc.

use "C:\Users\Freeze Factor\Documents\data\Democratic Accountability\countrylevel16.0.dta", clear

foreach y in dscldd nrd parpold sdd coposal_d {
forvalues x = 1/5 {
gsort -`y'`x'
gen `y'`x'_rank=_n
}
}


keep ccode rcode e1nwe e2nwe e3nwe e4nwe e5nwe b15_impnwe dwallnrnwe dwnpolwe dwnsdwe epgdppcppp eplitadult pcdemsat pcdemprefer pcethnof_relfrac pcethnof_ethfrac pctrust pcwvshappiness pecgeducexp pegini pehimartaxrt penews pepopdens pwomlegis
merge 1:1 ccode using "C:\Users\Freeze Factor\Documents\ConferencePapers\APSA2010Programmatic\apsa2010freezekitscheltcountrylevel.dta
drop if _merge==1
drop _merge

*Scatter of nrd1 and parpold1
twoway (scatter nrd1 parpold1 if rcode==1, mcolor(red) msymbol(circle) mlabel(ccodewb) mlabcolor(red)) (lfit nrd1 parpold1) (scatter nrd1 parpold1 if rcode==2, mcolor(black) msymbol(triangle) mlabel(ccodewb) mlabcolor(black)) (scatter nrd1 parpold1 if rcode==3, mcolor(forest_green) msymbol(lgx) mlabel(ccodewb) mlabcolor(forest_green)) (scatter nrd1 parpold1 if rcode==4, mcolor(orange) msymbol(circle_hollow) mlabel(ccodewb) mlabcolor(orange)) (scatter nrd1 parpold1 if rcode==5, mcolor(navy) msymbol(plus) mlabel(ccodewb) mlabcolor(navy)), ytitle(Redistribution Salience Score) xtitle(Redistribution Polarization Score) legend(on order(1 "Advanced Capitalist" 3 "Post-Communist" 4 "Latin America" 5 "Sub-Saharan Africa" 6 "Asia/Mideast"))

*Scatter of sdd1 and parpold1


*Scatter of b15_impnwe and coposal_4
twoway (scatter b15_impnwe coposal_4 if rcode==1, mcolor(red) msymbol(circle) mlabel(ccodewb) mlabcolor(red)) (lfit b15_impnwe coposal_4) (scatter b15_impnwe coposal_4 if rcode==2, mcolor(black) msymbol(triangle) mlabel(ccodewb) mlabcolor(black)) (scatter b15_impnwe coposal_4 if rcode==3, mcolor(forest_green) msymbol(lgx) mlabel(ccodewb) mlabcolor(forest_green)) (scatter b15_impnwe coposal_4 if rcode==4, mcolor(orange) msymbol(circle_hollow) mlabel(ccodewb) mlabcolor(orange)) (scatter b15_impnwe coposal_4 if rcode==5, mcolor(navy) msymbol(plus) mlabel(ccodewb) mlabcolor(navy)), ytitle(Additive Clientelism Index (B15)) xtitle(CoSalPo 4 issues) legend(on order(1 "Advanced Capitalist" 3 "Post-Communist" 4 "Latin America" 5 "Sub-Saharan Africa" 6 "Asia/Mideast"))

*Scatter of e2nwe and coposal_4
twoway (scatter coposal_4 e2nwe if rcode==1, mcolor(red) msymbol(circle) mlabel(ccodewb) mlabcolor(red)) (lfit coposal_4 e2nwe) (scatter coposal_4 e2nwe if rcode==2, mcolor(black) msymbol(triangle) mlabel(ccodewb) mlabcolor(black)) (scatter coposal_4 e2nwe if rcode==3, mcolor(forest_green) msymbol(lgx) mlabel(ccodewb) mlabcolor(forest_green)) (scatter coposal_4 e2nwe if rcode==4, mcolor(orange) msymbol(circle_hollow) mlabel(ccodewb) mlabcolor(orange)) (scatter coposal_4 e2nwe if rcode==5, mcolor(navy) msymbol(plus) mlabel(ccodewb) mlabcolor(navy)), ytitle(Operational Programmatic Score (CoSalPo 4)) xtitle(Symbolic Programmatic Evaluation(E2)) legend(on order(1 "Advanced Capitalist" 3 "Post-Communist" 4 "Latin America" 5 "Sub-Saharan Africa" 6 "Asia/Mideast"))

replace epgdppcppp=epgdppcppp/1000

*Scatter of GDP per capita and coposal_4
twoway (scatter epgdppcppp coposal_4 if rcode==1, mcolor(red) msymbol(circle) mlabel(ccodewb) mlabcolor(red)) (lfit epgdppcppp coposal_4) (scatter epgdppcppp coposal_4 if rcode==2, mcolor(black) msymbol(triangle) mlabel(ccodewb) mlabcolor(black)) (scatter epgdppcppp coposal_4 if rcode==3, mcolor(forest_green) msymbol(lgx) mlabel(ccodewb) mlabcolor(forest_green)) (scatter epgdppcppp coposal_4 if rcode==4, mcolor(orange) msymbol(circle_hollow) mlabel(ccodewb) mlabcolor(orange)) (scatter epgdppcppp coposal_4 if rcode==5, mcolor(navy) msymbol(plus) mlabel(ccodewb) mlabcolor(navy)), ytitle(GDPPC at PPP 1000s USD) xtitle(CoSalPo 4 issues) legend(on order(1 "Advanced Capitalist" 3 "Post-Communist" 4 "Latin America" 5 "Sub-Saharan Africa" 6 "Asia/Mideast"))

*Scatter of b15_impnwe and coposal_d1 (economic redistribution)
twoway (scatter b15_impnwe coposal_d1 if rcode==1, mcolor(red) msymbol(circle) mlabel(ccodewb) mlabcolor(red)) (lfit b15_impnwe coposal_d1) (scatter b15_impnwe coposal_d1 if rcode==2, mcolor(black) msymbol(triangle) mlabel(ccodewb) mlabcolor(black)) (scatter b15_impnwe coposal_d1 if rcode==3, mcolor(forest_green) msymbol(lgx) mlabel(ccodewb) mlabcolor(forest_green)) (scatter b15_impnwe coposal_d1 if rcode==4, mcolor(orange) msymbol(circle_hollow) mlabel(ccodewb) mlabcolor(orange)) (scatter b15_impnwe coposal_d1 if rcode==5, mcolor(navy) msymbol(plus) mlabel(ccodewb) mlabcolor(navy)), ytitle(Additive Clientelism Index) xtitle(Redistribution CoPolSal Score) legend(on order(1 "Advanced Capitalist" 3 "Post-Communist" 4 "Latin America" 5 "Sub-Saharan Africa" 6 "Asia/Mideast"))

*standardizing from 0 to 1

foreach x of varlist b15_impnwe-e5nwe  copolsal_mean- dissal_3econ {
egen `x'max=max(`x')
egen `x'min=min(`x')
gen `x'standard=`x'
replace `x'standard=`x'standard-`x'min
replace `x'standard=`x'standard/(`x'max-`x'min)
drop `x'max
drop `x'min
}


gen cps4e2diff=e2nwestandard-coposal_4standard

twoway (scatter coposal_4 e2nwestandard if rcode==1, mcolor(red) msymbol(circle) mlabel(ccodewb) mlabcolor(red)) (lfit coposal_4 e2nwestandard) (scatter coposal_4 e2nwestandard if rcode==2, mcolor(black) msymbol(triangle) mlabel(ccodewb) mlabcolor(black)) (scatter coposal_4 e2nwestandard if rcode==3, mcolor(forest_green) msymbol(lgx) mlabel(ccodewb) mlabcolor(forest_green)) (scatter coposal_4 e2nwestandard if rcode==4, mcolor(orange) msymbol(circle_hollow) mlabel(ccodewb) mlabcolor(orange)) (scatter coposal_4 e2nwestandard if rcode==5, mcolor(navy) msymbol(plus) mlabel(ccodewb) mlabcolor(navy)), ytitle(Operational Programmatic Score (CoSalPo 4)) xtitle(Symbolic Programmatic Evaluation(E2)) legend(on order(1 "Advanced Capitalist" 3 "Post-Communist" 4 "Latin America" 5 "Sub-Saharan Africa" 6 "Asia/Mideast"))

*Finding outliers - distance away from a 45 degree line in a scatterplot of two variables
local first e2nwestandard
local second coposal_4standard
gen diff=`first'-`second'
gen distance=(2*diff^2)^0.5

use "C:\Users\Freeze Factor\Documents\data\Democratic Accountability\partylevel19.0.dta", clear

keep ccode party e1 e2 e3 e4 e5 b11 b15 dw dwallnr dwsd rcode epgdppc

merge 1:1 ccode party using "C:\Users\Freeze Factor\Documents\ConferencePapers\APSA2010Programmatic\apsa2010freezekitscheltpartylevel.dta"

drop if _merge==1
drop _merge

pwcorr b15- e5  copolsal_mean- coposal_3same cosal_3econ cosal_3same cosal_4 cosal_8 cosal_mean d1 p31

*Scattergram between L-R position and Coposal 4
twoway (scatter coposal_4 dw if rcode==1, mcolor(red) msymbol(circle)) (lfit coposal_4 dw) (scatter coposal_4 dw if rcode==2, mcolor(black) msymbol(triangle)) (scatter coposal_4 dw if rcode==3, mcolor(forest_green) msymbol(lgx)) (scatter coposal_4 dw if rcode==4, mcolor(orange) msymbol(circle_hollow)) (scatter coposal_4 dw if rcode==5, mcolor(navy) msymbol(plus)), ytitle(Coposal 4) xtitle(Left-Right position) legend(on order(1 "Advanced Capitalist" 3 "Post-Communist" 4 "Latin America" 5 "Sub-Saharan Africa" 6 "Asia/Mideast"))

*standardizing from 0 to 1

foreach x of varlist b15- e5  copolsal_mean- coposal_3same cosal_3econ cosal_3same cosal_4 cosal_8 cosal_mean {
egen `x'max=max(`x')
egen `x'min=min(`x')
gen `x'standard=`x'
replace `x'standard=`x'standard-`x'min
replace `x'standard=`x'standard/(`x'max-`x'min)
drop `x'max
drop `x'min
}

gen cps4e2diff=e2-coposal_4

reg cps4e2diff b15 epgdppc p31

reg cps4e2diff b15 epgdppc p31 dwsd

reg cps4e2diff b15 epgdppc p31 dwsd coposal_3econ

*Finding outliers - distance away from a 45 degree line in a scatterplot of two variables
local first e2standard
local second copolsal_meanstandard
gen diff=`first'-`second'
gen distance=(2*diff^2)^0.5

*/


