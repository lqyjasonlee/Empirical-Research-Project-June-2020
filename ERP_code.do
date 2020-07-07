# delimit ;
set more off;

clear all ;
set maxvar 20000 ;
set matsize 11000 ;

capture log close ;
log using erp_result_new.txt , replace text ;

use ERP_data_qiyuan.dta ;


* Normalized trend variable ;
gen trend = (year - 1980)/12 ;

tsset statenum year ;

* Estimate baseline model in first-differences ;
xi: reg D.lpc_viol D.irca D.lpop D.UR D.povrate D.ATWI D.crack_index D.la14pc D.opc D.logincome D.logemp i.year , cluster(statenum) ;
xi: reg D.lpc_prop D.irca D.lpop D.UR D.povrate D.ATWI D.crack_index D.la14pc D.opc D.logincome D.logemp i.year , cluster(statenum) ;
xi: reg D.lpc_all D.irca D.lpop D.UR D.povrate D.ATWI D.crack_index D.la14pc D.opc D.logincome D.logemp i.year , cluster(statenum) ;

* Generate variables for LASSO ;

local tdums = "_Iyear_1982 _Iyear_1983 _Iyear_1984 _Iyear_1985 _Iyear_1986 _Iyear_1987 _Iyear_1988 _Iyear_1989 _Iyear_1990 _Iyear_1991 _Iyear_1992 _Iyear_1993 _Iyear_1994 _Iyear_1995 _Iyear_1996 _Iyear_1997 _Iyear_1998 _Iyear_1999" ;
local xx = "lpop UR povrate ATWI crack_index la14pc opc logincome logemp" ;

* Differences ;
local Dxx ;
foreach x of local xx { ;
	gen D`x' = D.`x' ;
	local tempname = "D`x'" ;
	local Dxx : list Dxx | tempname ;
} ;

* Squared Differences
local Dxx2 ;
foreach x of local Dxx { ;
	gen `x'2 = `x'^2 ;
	local tempname = "`x'2" ;
	local Dxx2 : list Dxx2 | tempname ;
} ;

* Difference Interactions
local DxxInt ;
local nxx : word count `Dxx' ;
forvalues ii = 1/`nxx' { ;
	local start = `ii'+1 ;
	forvalues jj = `start'/`nxx' { ;
		local temp1 : word `ii' of `Dxx' ;
		local temp2 : word `jj' of `Dxx' ;
		gen `temp1'X`temp2' = `temp1'*`temp2' ;
		local tempname = "`temp1'X`temp2'" ;
		local DxxInt : list DxxInt | tempname ;
	} ;
} ;		

* Lags ;
local Lxx ;
foreach x of local xx { ;
	gen L`x' = L.`x' ;
	local tempname = "L`x'" ;
	local Lxx : list Lxx | tempname ;
} ;

* Squared Lags
local Lxx2 ;
foreach x of local Lxx { ;
	gen `x'2 = `x'^2 ;
	local tempname = "`x'2" ;
	local Lxx2 : list Lxx2 | tempname ;
} ;

* Means ;
local Mxx ;
foreach x of local xx { ;
	by statenum: egen M`x' = mean(`x') ;
	local tempname = "M`x'" ;
	local Mxx : list Mxx | tempname ;
} ;

* Squared Means ;
local Mxx2 ;
foreach x of local Mxx { ;
	gen `x'2 = `x'^2 ;
	local tempname = "`x'2" ;
	local Mxx2 : list Mxx2 | tempname ;
} ;

* Initial Levels ;
local xx0 ;
foreach x of local xx { ;
	by statenum: gen `x'0 = `x'[1] ;
	local tempname = "`x'0" ;
	local xx0 : list xx0 | tempname ;
} ;

* Squared Initial Levels ;
local xx02 ;
foreach x of local xx0 { ;
	gen `x'2 = `x'^2 ;
	local tempname = "`x'2" ;
	local xx02 : list xx02 | tempname ;
} ;

* Initial Differences ;
local Dxx0 ;
foreach x of local Dxx { ;
	by statenum: gen `x'0 = `x'[2] ;
	local tempname = "`x'0" ;
	local xx0 : list xx0 | tempname ;
} ;

* Squared Initial Differences ;
local Dxx02 ;
foreach x of local Dxx0 { ;
	gen `x'2 = `x'^2 ;
	local tempname = "`x'2" ;
	local Dxx02 : list Dxx02 | tempname ;
} ;

* Interactions with trends ;
local biglist : list Dxx | Dxx2 ;
local biglist : list biglist | DxxInt ;
local biglist : list biglist | Lxx ;
local biglist : list biglist | Lxx2 ;
local biglist : list biglist | Mxx ;
local biglist : list biglist | Mxx2 ;
local biglist : list biglist | xx0 ;
local biglist : list biglist | xx02 ;
local biglist : list biglist | Dxx0 ;
local biglist : list biglist | Dxx02 ;

local IntT ;
local nxx : word count `biglist' ;
foreach x of local biglist { ;
	gen `x'Xt = `x'*trend ;
	gen `x'Xt2 = `x'*(trend^2) ;
	local tempname = "`x'Xt `x'Xt2" ;
	local IntT : list IntT | tempname ;
} ;
	
local shared : list biglist | IntT ;

* Violence specific controls ;
gen Dirca = D.irca ;
by statenum: gen viol0 = irca[1] ;
by statenum: gen Dviol0 = Dirca[2] ;
gen viol02 = viol0^2 ;
gen Dviol02 = Dviol0^2 ;
gen viol0Xt = viol0*trend ;
gen viol0Xt2 = viol0*(trend^2) ;
gen viol02Xt = viol02*trend ;
gen viol02Xt2 = viol02*(trend^2) ;
gen Dviol0Xt = Dviol0*trend ;
gen Dviol0Xt2 = Dviol0*(trend^2) ;
gen Dviol02Xt = Dviol02*trend ;
gen Dviol02Xt2 = Dviol02*(trend^2) ;

local contviol = "viol0 viol0Xt viol0Xt2 viol02 viol02Xt viol02Xt2 
			Dviol0 Dviol0Xt Dviol0Xt2 Dviol02 Dviol02Xt Dviol02Xt2" ;

local AllViol : list contviol | shared ;
			
* Property specifc controls ;
* gen Dprop = D.efaprop ;
by statenum: gen prop0 = irca[1] ;
by statenum: gen Dprop0 = Dirca[2] ;
gen prop02 = prop0^2 ;
gen Dprop02 = Dprop0^2 ;
gen prop0Xt = prop0*trend ;
gen prop0Xt2 = prop0*(trend^2) ;
gen prop02Xt = prop02*trend ;
gen prop02Xt2 = prop02*(trend^2) ;
gen Dprop0Xt = Dprop0*trend ;
gen Dprop0Xt2 = Dprop0*(trend^2) ;
gen Dprop02Xt = Dprop02*trend ;
gen Dprop02Xt2 = Dprop02*(trend^2) ;

local contprop = "prop0 prop0Xt prop0Xt2 prop02 prop02Xt prop02Xt2 
			Dprop0 Dprop0Xt Dprop0Xt2 Dprop02 Dprop02Xt Dprop02Xt2" ;
			
local AllProp : list contprop | shared ;

* All crime specific controls ;
gen Dall = D.irca ;
by statenum: gen all0 = irca[1] ;
by statenum: gen Dall0 = Dall[2] ;
gen all02 = all0^2 ;
gen Dall02 = Dall0^2 ;
gen all0Xt = all0*trend ;
gen all0Xt2 = all0*(trend^2) ;
gen all02Xt = all02*trend ;
gen all02Xt2 = all02*(trend^2) ;
gen Dall0Xt = Dall0*trend ;
gen Dall0Xt2 = Dall0*(trend^2) ;
gen Dall02Xt = Dall02*trend ;
gen Dall02Xt2 = Dall02*(trend^2) ;

local contall = "all0 all0Xt all0Xt2 all02 all02Xt all02Xt2 
			Dall0 Dall0Xt Dall0Xt2 Dall02 Dall02Xt Dall02Xt2" ;
			
local AllAll : list contall | shared ;

* Differenced outcomes ;
gen Dyviol = D.lpc_viol ;
gen Dyprop = D.lpc_prop ;
gen Dyall = D.lpc_all ;
			
drop if trend == 0 ;

* Regression using naive method ;
reg Dyviol irca `AllViol' `tdums' , cluster(statenum) ;
reg Dyprop irca `AllProp' `tdums' , cluster(statenum) ;
reg Dyall irca `AllAll' `tdums' , cluster(statenum) ;

* Using selection method ;

* Violence Outcome ;
lassoShooting Dyviol `AllViol' , controls(`tdums') lasiter(100) verbose(0) fdisplay(0) ;
local yvSel `r(selected)' ;
di "`yvSel'" ;

* Violence IRCA ;
lassoShooting irca `AllViol' , controls(`tdums') lasiter(100) verbose(0) fdisplay(0) ;
local xvSel `r(selected)' ;
di "`xvSel'" ;

* Get union of selected instruments ;
local vDS : list yvSel | xvSel ;

* Violence equation with selected controls ;
reg Dyviol irca `vDS' `tdums' , cluster(statenum) ;


* Property Outcome ;
lassoShooting Dyprop `AllProp' , controls(`tdums') lasiter(100) verbose(0) fdisplay(0) ;
local ypSel `r(selected)' ;
di "`ypSel'" ;

* Property IRCA ;
lassoShooting irca `AllProp' , controls(`tdums') lasiter(100) verbose(0) fdisplay(0) ;
local xpSel `r(selected)' ;
di "`xpSel'" ;

* Get union of selected instruments ;
local pDS : list ypSel | xpSel ;

* Property equation with selected controls ;
reg Dyprop irca `pDS' `tdums' , cluster(statenum) ;


* All Outcome ;
lassoShooting Dyall `AllAll' , controls(`tdums') lasiter(100) verbose(0) fdisplay(0) ;
local yaSel `r(selected)' ;
di "`yaSel'" ;

* All IRCA ;
lassoShooting irca `AllAll' , controls(`tdums') lasiter(100) verbose(0) fdisplay(0) ;
local xaSel `r(selected)' ;
di "`xaSel'" ;

* Get union of selected instruments ;
local aDS : list yaSel | xaSel ;

* Property equation with selected controls ;
reg Dyall irca `aDS' `tdums' , cluster(statenum) ;


log close ;

