
set more off
capture program drop datagen
global cloudstor "C:\Users\z3312911\Cloudstor"
sysdir set PLUS "$cloudstor\ado\"

program datagen

scalar drop _all

// Define scalars which are passed to the program by the simulation code
scalar misp="`1'" // Type of misspecification (if any)
local numdata="`2'" // Number of simulation iterations

if misp=="None" {
	global scenario "1"
}
if misp=="Exposure" {
	global scenario "2"
}
if misp=="Outcome" {
	global scenario "3"
}
if misp=="Both" {
	global scenario "4"
}

	
forvalues dat=1/`numdata' {
	clear
	set obs 250
	gen id=_n

	matrix C1=(1,0.3,0.3 \ 0.3,1,0.3 \ 0.3,0.3,1)
	matrix C2=(1,0.9,0.8,0.7,0.6 \ 0.9,1,0.9,0.8,0.7 \ 0.8,0.9,1,0.9,0.8 \ 0.7,0.8,0.9,1,0.9 \ 0.6,0.7,0.8,0.9,1)

	matrix othco=(2.0,1.5,0.5,-0.8,0.5)
	matrix parco=(1.5,1.0,0.0,0.5,-0.5,0.7,0.3)
	matrix outco=(0.5,1.0,1.0,0.0,1.0,0.5,0.8,-0.5)

	drawnorm oa ob oc, corr(C1) means(0 0 0) sds(1 1 1)
	drawnorm u0 u1 u2 u3 u4, corr(C2) means(0 0 0 0 0) sds(1 1 1 1 1)

	gen el_ijk=rnormal(0,1)
	gen ea_ijk=rnormal(0,1)
	gen ey_ijk=rnormal(0,1)
	
	reshape long u, i(id) j(obs)
	xtset id obs
	
	gen l=.
	gen a=.
	gen y=.
	gen la=rbinomial(1,0.2)
	gen ll=rbinomial(1,0.25)

	forvalues i=0/4 {

		replace la=l1.a if obs==`i' & l1.a!=.
		replace ll=l1.l if obs==`i' & l1.l!=.
		
		replace l = runiform()<invlogit(-3.7 + ///
		othco[1,1]*la + othco[1,2]*ll + ///
		othco[1,3]*oa + othco[1,4]*ob + othco[1,5]*oc + ///
		3*u + el_ijk) if obs==`i'

			if misp=="Exposure" | misp=="Both" { // Interaction type misspec of propensity model with interaction defined a_j*l_i
				matrix parco[1,3]=0.5
			}
			
		replace a = runiform()<invlogit(-2.3 + ///
		parco[1,1]*la + parco[1,2]*l + parco[1,3]*obs + parco[1,4]*ll + ///
		parco[1,5]*oa + parco[1,6]*ob + parco[1,7]*oc + ea_ijk) if obs==`i'

			if misp=="Outcome" | misp=="Both" { // Interaction type misspec of propensity model with interaction defined a_j*l_i
				matrix outco[1,4]=-0.5
			}
			
		replace y = rnormal(5 + ///
		outco[1,1]*la + outco[1,2]*a + outco[1,3]*l + outco[1,4]*obs + ///
		outco[1,5]*ll + outco[1,6]*oa + outco[1,7]*ob + outco[1,8]*oc + ///
		3*u + ey_ijk,1) if obs==`i'
		
	}

	drop if obs==0
	keep id y a l la ll oa ob oc obs
	order obs, after(id)
	order y, after(obs)
	order a, after(y)
	order la, after(a)
	order l, after(la)
	order ll, after(l)
	
	save "$cloudstor\PhD\Paper 3 - Methods simulation\Alternate Simulations\Scenario B$scenario\Simulation $scenario-`dat'.dta", replace

}
		
end

local numdata=500

set seed 349447
datagen None `numdata'
set seed 349447
datagen Exposure `numdata'
set seed 349447
datagen Outcome `numdata'
set seed 349447
datagen Both `numdata'

/*
capture program drop getscalars
program getscalars

global i=$i+1
use "$cloudstor\PhD\Paper 3 - Methods simulation\Alternate Simulations\Scenario B1\Simulation 1-$i.dta", clear
su a
	scalar lprop1=r(mean)
use "$cloudstor\PhD\Paper 3 - Methods simulation\Alternate Simulations\Scenario B2\Simulation 2-$i.dta", clear
su a
	scalar lprop2=r(mean)
use "$cloudstor\PhD\Paper 3 - Methods simulation\Alternate Simulations\Scenario B3\Simulation 3-$i.dta", clear
su a
	scalar lprop3=r(mean)
use "$cloudstor\PhD\Paper 3 - Methods simulation\Alternate Simulations\Scenario B4\Simulation 4-$i.dta", clear
su a
	scalar lprop4=r(mean)
end

global i=0
simulate lprop1=lprop1 lprop2=lprop2 lprop3=lprop3 lprop4=lprop4, reps(500): getscalars
su

