*! rfbunch version date 20211130
* Author: Martin Eckhoff Andresen
* This program is part of the rfbunch package.
cap prog drop rfbunch
program rfbunch, eclass sortpreserve
	syntax varlist(min=1) [if] [in],  CUToff(real) bw(real) [ ///
	LIMits(numlist min=1 max=2 >=0) ///
	notch(numlist min=2 max=3 >=0) ///
	kink(numlist min=2 max=2 >0) ///
	tcr(numlist min=3 max=3) ///
	init(numlist min=1 max=2 >0) ///
	POLynomial(numlist min=0 integer >=0) ///
	localbw(string) ///
	localkernel(string) ///
	ADJust(string) ///
	nofill ///
	constant ///
	placebo ///
	xtype(numlist <=3 >=0) ///
	]
	
	quietly {

		//check options
		gettoken varlist xvars: varlist
		
		tempname polynomials
		loc numxvars: word count `xvars'
		loc numpoly: word count `polynomial'
		
		if `numpoly'>`numxvars'+1 {
			noi di as error "Specify no more numbers than the number of variables in polynomial()."
			exit 301
		}
		if `numpoly'==0 {
			mat `polynomials'=7
			foreach xvar in `xvars' {
				mat `polynomials'=`polynomials' \ 1
				}
			}	
		else {
			foreach pol in `polynomial' {
				mat `polynomials'=nullmat(`polynomials') \ `pol'
				loc lastpol=`pol'
				}
				if `numpoly'<`numxvars'+1&`numpoly'>1 {
					forvalues k=1/`=`numxvars'+1-`=rowsof(`polynomials')'' {
						mat `polynomials'=`polynomials' \ `lastpol'
					}
				}
			}
		
		mat rownames `polynomials'=`varlist' `xvars'
		mat colnames `polynomials'=polynomial
		
		if "`xvars'"!="" {
			tempvar xtypes
			if "`xtype'"=="" {
				foreach var in `xvars' {
					mat `xtypes'=nullmat(`xtypes') \ 0
					}
				}
			else {
				local numxvals: word count `xtype'
				if `numxvals'>`numxvars' {
					noi di as error "Specify no more numbers than the number of x variables in numlist xtypes()."
					exit 301
				}
				foreach xval in `xtype' {
					mat `xtypes'=nullmat(`xtypes') \ `xval'
					loc lastxval=`xval'
				}
				
				forvalues k=1/`=`numxvars'-`numxvals'' {
					mat `xtypes'=nullmat(`xtypes') \ `lastxval'
				}
			}
		
		mat rownames `xtypes'=`xvars'
		mat colnames `xtypes'=xtypes
		}

		

		
		cap which moremata.hlp
		
		if _rc!=0 {
			noi di in red "moremata needed. Install using "ssc install moremata"".
			exit 301
			}
		
		if "`adjust'"!="" {
			if !inlist("`adjust'","y","x","none","logx") {
				noi di in red "Option adjust can only take values y (Chetty et al.), x and logx (Andresen and Thorvaldsen) or none"
				exit 301
				}
			}
		if "`adjust'"=="x" loc type=2
		else if "`adjust'"=="logx" loc type=3
		else if "`adjust'"=="y" loc type=1
		else loc type=0
		
		if "`limits'"!="" {
			gettoken L H: limits
			if "`H'"=="" loc H=0
		}
		else {
			loc L=1
			loc H=0
		}
		
		loc zL=`cutoff'-`L'*`bw'
		loc zH=`cutoff'+`H'*`bw'
		
		if "`fill'"=="nofill" loc fill=0
		else loc fill=1
		
		if "`localbw'"!="" {
			if "`localbw'"!="rdbw" {
				cap confirm scalar `localbw'
				if _rc!=0 {
					noi di in red "Option localbw() can only take the values rdbw or a single nonnegative number."
					exit 301
				}
				else {
					if `localbw'<0 {
						noi di in red "Option localbw() cannot be negative."
						exit 301
					}
				}
			}
			else {
				cap which rdbwselect
				if _rc!=0 {
						noi di in red "Option localbw(rdbw) requires the rdrobust package."
						exit 301
				}
			}
		}
		
		if !inlist("`localkernel'","","triangular","epanechnikov","uniform") {
				noi di in red "Option localkernel() can only be triangular, epanechnikov or uniform."
				exit 301
		}
		
		//check tcr, kink, notch options
		if "`tcr'"!="" {
			if "`adjust'"!="x" {
				noi di as text "Option tcr() requires adjust(x) for identification. Option tcr ignored."
				loc tcr
			}
			else {
				gettoken t_tcr tcr: tcr
				gettoken eta_tcr r_tcr: tcr
				if !inrange(`t_tcr',0,1)|!inrange(`eta_tcr',0,1) {
					noi di as error "syntax of tcr is tcr(t eta r), where t and eta are numbers between 0 and 1"
					exit 301
				}
				if "`init'"=="" {
					tempname initmat
					mat `initmat'=(0.2,0.2)
				}
				else {
					loc numinit: word count `init'
					if `numinit'!=2 {
						noi di as error "Specify two number in init() when using tcr() - mu,epsilon"
						exit 301
					}
					gettoken mu epsilon: init
					tempname initmat
					mat `initmat'=(`mu',`epsilon')
					}
				}
			}
	
		if "`notch'"!="" {
			gettoken t0_notch rest: notch
			gettoken t1_notch rest: rest
			if "`rest'"=="" loc deltaT_notch=0
			else loc deltaT_notch=`rest'
			if !inrange(`t0_notch',0,1)|!inrange(`t1_notch',0,1)|`deltaT_notch'<0 {
				noi di in red "syntax of notch is notch(t0 t1 [deltaT]), where t0 and t1 are tax rates between 0 and 1 and deltaT are nonnegative."
				exit 301
			}
			if "`init'"=="" loc init=0.1
			else {
				loc numinit: word count `init'
				if `numinit'!=1 {
					noi di as error "Specify only one number in init() when using notch()"
					exit
				}
			}
		}
		
		if "`kink'"!="" {
			gettoken t0_kink t1_kink: kink
			if !inrange(`t0_kink',0,1)|!inrange(`t1_kink',0,1) {
				noi di in red "syntax of kink is kink(t0 t1), where t0 and t1 are tax rates between 0 and 1."
				exit 301
			}
			if "`init'"=="" loc init=0
			else {
				loc numinit: word count `init'
				if `numinit'!=1 {
					noi di as error "Specify only one number in init() when using kink()"
					exit
				}
			}
		}
		
		count if `varlist'<`cutoff'
		if r(N)==0 {
			noi di as error "No individuals in sample allocates below cutoff."
			exit 301
			}
		
		count if `varlist'>`cutoff'
		if r(N)==0 {
			noi di as error "No individuals in sample allocates above cutoff."
			exit 301
			}
		
		tempvar resid freq0 freq touse useobs bin
		tempname table cutvals cf b adj_table obsbins shiftmat
		marksample touse
		preserve
		drop if !`touse'
		keep `varlist'
		
		loc N=_N
		count if `varlist'>`zL'
		loc BM=r(N)
		count if `varlist'>`zL'&`varlist'<=`zH'
		loc Bunchmass=r(N)
		
		loc colfreq `varlist' frequency
		gen `useobs' = `varlist'<=`zL'|`varlist'>`zH'
		
		su `varlist' if `varlist'>`cutoff'
		if ((r(min)>`cutoff'+`bw')) {
			noi di as text "Note: Hole detected above cutoff. Bins adjusted above cutoff to avoid half-empty bins."
			loc hole=1
			loc minabove=r(min)
		} 
		else loc hole=0
		if "`adjust'"!=""|`hole'==1 {
			loc minabove=r(min)
		}
		else loc minabove=`cutoff'
		
		forvalues i=1/`=`polynomials'[1,1]' {
			if "`rhsvars'"=="" loc rhsvars  c.`varlist'
			else loc rhsvars `rhsvars'##c.`varlist'
			mat `cutvals'=nullmat(`cutvals') \ `cutoff'^`i'
			loc coleq `coleq' counterfactual_frequency
		}
	
		loc coleq `coleq' counterfactual_frequency
		mat `cutvals'=1 \ `cutvals'
		
		mata: st_matrix("`table'",fill(st_data(.,"`varlist'"),`bw',`cutoff',`zH',0,`type',0,`cutoff',`hole'))
		//Get counterfactual and adjust, if using
		if inlist("`adjust'","x","y","logx") {
			mata: shift=shifteval(st_data(selectindex(st_data(.,"`useobs'")),"`varlist'"),`zL',`zH',`=`polynomials'[1,1]',`BM',`bw',`type',10,0,`fill',`cutoff',`hole')
			mata: st_matrix("`b'",shift)
			mata: st_matrix("`adj_table'",fill(st_data(.,"`varlist'"),`bw',`cutoff',`cutoff',`=`b'[1,`=colsof(`b')']',`type',0,`cutoff',`hole'))
			mat `cf'=`b'[1,1..`=colsof(`b')-1']
			mat `b'=`b'[1,2..`=colsof(`b')-1'],`b'[1,1],`b'[1,`=colsof(`b')']
			local shift=`b'[1,`=colsof(`b')']
			fvexpand `rhsvars'
			loc names `r(varlist)' _cons shift
			loc coleq `coleq' bunching
			loc adjnames adj_bin adj_freq
			}
		else {
			mata: data=fill(st_data(selectindex(st_data(.,"`useobs'")),"`varlist'"),`bw',`zL',`zH',0,0,`fill',`cutoff',`hole')
			mata: xbin=J(rows(data[.,1]),1,1)
			mata: for (p=1; p<=`=`polynomials'[1,1]'; p++) xbin=xbin,data[.,1]:^p
			mata: b=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
			mata: st_matrix("`b'",b)
			mat `cf'=`b'
			mat `b'=`b'[1,2..`=colsof(`b')'],`b'[1,1]
			fvexpand `rhsvars'
			loc names `r(varlist)' _cons
			local shift=0
		}
		
		tempname freq0b freq0tau
		mata: st_matrix("`freq0b'",(polyeval(polyinteg(st_matrix("`cf'"),1),`zH')-polyeval(polyinteg(st_matrix("`cf'"),1),`zL'))/`bw')
		mata: st_matrix("`freq0tau'",(polyeval(st_matrix("`cf'"),`cutoff')))
		loc B=`=`Bunchmass'-`freq0b'[1,1]'
		mat `b'=`b',`B',`=`B'/`N'',`=`B'/`freq0b'[1,1]',`=`B'/`freq0tau'[1,1]' //Number of bunchers, bunchers share of sample, normalized bunching, excess mass
		loc coleq `coleq' bunching bunching bunching bunching
		loc names `names' number_bunchers share_sample normalized_bunching excess_mass 
		
		if "`constant'"!="constant" {
			mata: meannonbunch=(polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff') -polyeval(polyinteg((0,st_matrix("`cf'")),1),`zL'))/(polyeval(polyinteg((st_matrix("`cf'")),1),`cutoff') -polyeval(polyinteg((st_matrix("`cf'")),1),`zL'))
			mata: st_numscalar("meannonbunch",meannonbunch)
			}
		else scalar meannonbunch=(`cutoff'-`zL')/2
		
		if `B'<0|"`placebo'"!="" {
			if `B'<0&"`placebo'"=="" noi di as text "Negative estimates of B - no bunching in the bunching region. Marginal response, total response and counterfactual mean among bunchers cannot be calculated."
			mat `b'=`b',meannonbunch
			loc coleq `coleq' bunching
			loc names `names' mean_nonbunchers
			
			}
			
		else {
			if "`constant'"!="constant" {
				mata: eresp=eresp(`B',`cutoff',st_matrix("`cf'"),`bw')
				if eresp==. {
					noi di as error 	"No real roots found above the cutoff to the integral of the counterfactual."
					noi di as error 	"Marginal response and various other parameters cannot be estimated. You could"
					noi di as error		"try using the constant approximation to the density in the bunching region by"
					noi di as error		"specifying the option constant."
					exit 301
				}
				mata: st_numscalar("eresp",eresp)
				mata: meanbunch=(polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff'+eresp) -polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff'))/(`bw'*`B')-(`cutoff')
				mata: st_numscalar("meanbunch",meanbunch)
				mata: totalresponse=(1/`bw')*(polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff'+eresp) -polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff'))-(`cutoff'*`B')
				mata: st_numscalar("totalresponse",totalresponse)
				}
			else {
				tempname predcut
				mat `predcut'=`cf'*`cutvals'
				scalar eresp=(`bw'*`B')/`predcut'[1,1]
				scalar meanbunch=eresp*`predcut'[1,1]-`B'*`cutoff'/2
				scalar totalresponse=eresp*`predcut'[1,1]
				}
			
			loc coleq `coleq' bunching bunching bunching bunching
			loc names `names' mean_nonbunchers marginal_response total_response average_response
			mat `b'=`b',meannonbunch,eresp,totalresponse,meanbunch
			}

		
		
		restore
	
		//ESTIMATE REPONSE ALONG OTHER ENDOGENOUS VARS
		if "`xvars'"!="" {
			if `B'<0&"`placebo'"=="" noi di as text "Mean counterfactual for bunchers and difference between bunchers means and this quantity cannot be calculated for alternative because B<0."
			
			preserve
			drop if !`touse'
			keep `varlist' `xvars'
			gen `useobs'=`varlist'<=`zL'|`varlist'>`zH'
			count if !`useobs'
			loc N_bunchreg=r(N)
			
			if "`localbw'"!="" {
				if "`localbw'"!="rdbw" {
					if `localbw'==0 {
						su `varlist'
						loc bwlow=`cutoff'-r(min)
						loc bwhi=r(max)-`cutoff'
					}
					else {
						loc bwlow=`localbw'
						loc bwhi=`localbw'
					}
				
					if inlist("`localkernel'","","triangular") {
						gen double w=1-abs(`varlist'-`cutoff')/`bwlow' if `varlist'>`cutoff'-`bwlow'&`varlist'<=`cutoff'
						replace w=1-abs(`varlist'-`cutoff')/`bwhi' if `varlist'<`cutoff'+`bwhi'&`varlist'>`cutoff'
					}
					else if "`localkernel'"=="epanechnikov" {
						gen double w=(3/4)*max((1-((`cutoff'-`varlist')/`bwlow')^2),0) if `varlist'>`cutoff'-`bwlow'&`varlist'<=`cutoff'
						replace w=(3/4)*max((1-((`cutoff'-`varlist')/`bwhi')^2),0) if `varlist'<`cutoff'+`bwhi'&`varlist'>`cutoff'
					}
					else {
						gen double w=1/2 `varlist'>`cutoff'-`bwlow'&`varlist'<=`cutoff'
						replace w=1/2 if `varlist'<`cutoff'+`bwhi'&`varlist'>`cutoff'
						}
					}
				loc localweights [aw=w]
				}
									
			tempname means integerbin adjustbin adjustz
			tempvar predy f0 f1 f above fabove mean_b_cf
			
			if `type'<2&`hole'==0 gen `bin'=ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2
			else gen `bin'=(`varlist'<=`cutoff')*(ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2) + (`varlist'>`cutoff')*(floor((`varlist'-`minabove'+2^-23)/`bw')*`bw'+`minabove'+`bw'/2)			
			sort `bin'
			egen `integerbin'=group(`bin')
			
			if `type'==2 {
				replace `bin'=(`varlist'<=`cutoff')*(ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2) + ((`varlist'>`cutoff')*(floor(((`varlist'-`minabove'+2^-23)/(1+`shift'))/`bw')*`bw'+`minabove'/(1+`shift')+`bw'/2))
				gen `adjustz'=`varlist'/(1+`shift'*(`varlist'>`cutoff'))
				}
			else if `type'==3 {
				replace `bin'=(`varlist'<=`cutoff')*(ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2) + ((`varlist'>`cutoff')*(floor(-log(1+`shift')+(`varlist'-`minabove'+2^-23)/`bw')*`bw'+`minabove'-log(1+`shift')+`bw'/2))
				gen `adjustz'=`varlist'-log(1+`shift')*(`varlist'>`cutoff')				
				}
			egen `adjustbin'=group(`bin')
			
			gen `above'=`varlist'>`cutoff'
			loc i=0
			foreach var in `xvars' {
				loc ++i
				
				if "`localbw'"=="rdbw" {
					rdbwselect `var' `varlist' if `useobs', c(`cutoff') kernel(`localkernel')
					cap drop w
					gen w=1-abs(`varlist'-`cutoff')/e(h_mserd) if `useobs' & inrange(`varlist',`cutoff'-e(h_mserd),`cutoff'+e(h_mserd))
					loc localweights [aw=w]
					mat `localbws'=nullmat(`localbws') \e(h_mserd)
					
					}
					
				if `=`polynomials'[`=`i'+1',1]'>0 {
					loc rhsvarsnl
					loc rhsvars
					forvalues k=1/`=`polynomials'[`=`i'+1',1]' {
						if `xtypes'[`i',1]==0|inlist("`adjust'","none","","y") loc xvar `varlist'
						else loc xvar `adjustz'
						if `k'==1 loc rhsvars c.`xvar'
						else loc rhsvars `rhsvars'##c.`xvar'
						if `xtypes'[`i',1]==1 {
								loc rhsvarsnl `rhsvarsnl' {beta`k'}*`xvar'^`k'*(1+{gamma}*`above')+
						}
					}
				}
				
				if `xtypes'[`i',1]==0  reg `var' `rhsvars' 1.`above' 1.`above'#(`rhsvars')  `localweights'  if `useobs' //no link
				else if `xtypes'[`i',1]==1  {
					reg `var' `rhsvars' if `varlist'<`zL'
					tempname init
					mat `init'=e(b)
					loc initials
					forvalues k=1/`=colsof(`init')-1' {
						loc initials `initials' beta`k' `=`init'[1,`k']'
					}
					nl (`var'= `rhsvarsnl' {beta0}*(1+{gamma}*`above') ) if `useobs', initial(beta0 `=_b[_cons]' `initials' gamma 0)
					fvexpand `rhsvars'
					loc names `names' `=subinstr("`=subinstr("`r(varlist)'","`adjustz'","`varlist'",.)'","1.`above'","above",.)' above _cons
				}
				else if `xtypes'[`i',1]==2 reg `var' `rhsvars' 1.`above' `localweights'  if `useobs' //assume x0/x1=constant, y in logs
				else if `xtypes'[`i',1]==3 reg `var' `rhsvars'  `localweights' if `useobs' //characterize, assume x0==x1
				
				mat `b'=`b',e(b)
				mat `f'=e(b)
				
				if `=`polynomials'[`=`i'+1',1]'>0 mat `f'=`f'[1,`=colsof(`f')'],`f'[1,1..`=`polynomials'[`=`i'+1',1]']
				else mat `f'=_b[_cons]

				mata:mean_nonbunchers=(polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`cutoff')-polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`zL'))/(polyeval(polyinteg(st_matrix("`cf'"),1),`cutoff')-polyeval(polyinteg(st_matrix("`cf'"),1),`zL'))	
				mata: st_numscalar("mean_nonbunchers",mean_nonbunchers)
				local colnames: colnames e(b)
				if `xtypes'[`i',1]!=1 loc names `names' `=subinstr("`=subinstr("`colnames'","`adjustz'","`varlist'",.)'","1.`above'","above",.)'
				forvalues j=1/`=colsof(e(b))' {
					loc coleq `coleq' `var'
				}
				predict double `predy' if inrange(`varlist',`zL',`cutoff'), residuals
				su `predy'
				loc excess=r(sum)
				drop `predy'
				mat `b'=`b',`excess',mean_nonbunchers,`excess'/`N_bunchreg'
				
				loc names `names' excess_value mean_nonbunchers average_prediction_error
				loc coleq `coleq' `var'_means `var'_means `var'_means
				if `B'>0&"`placebo'"=="" {
					mata: mean_counterfactual=(polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`cutoff'+`=eresp')-polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`cutoff'))/(polyeval(polyinteg(st_matrix("`cf'"),1),`cutoff'+`=eresp')-polyeval(polyinteg(st_matrix("`cf'"),1),`cutoff'))
					mata: st_numscalar("mean_b_cf",mean_counterfactual)
					mat `b'=`b',`=`excess'/`B'+mean_nonbunchers',mean_b_cf,`=`excess'/`B'+mean_nonbunchers-mean_b_cf'
					loc names `names' mean_bunchers mean_bunchers_cf bunchers_diff
					loc coleq `coleq' `var'_means `var'_means `var'_means 
					}
					
				reg `var' ibn.`integerbin', nocons
				mat `means'=e(b)
				
				mat `table'=`table',`means''
				loc colfreq `colfreq' `var'
								
				if `type'>=2 {
					reg `var' ibn.`adjustbin', nocons
					mat `means'=e(b)
					
					mat `adj_table'=`adj_table',`means''
					loc adjnames `adjnames' adj_`var'
				}
				
				
				}
			restore
			}
		
		//ESTIMATE ELASTICITIES
		tempname e
		if "`tcr'"!="" {
			mata: e=tcr(`t_tcr',`cutoff',`eta_tcr',`r_tcr',`=eresp',`shift',st_matrix("`initmat'"))
			mata: st_matrix("`e'",e)
			if `e'[1,5]!=0 {
				noi di as error "Error code `=errorcode' during numeric optimization using tcr(), se help mata optimize"
				noi di as error "Elasticity estimates not reported"
				mat `b'=`b',.,.,.,. //elasticity
			}
			else mat `b'=`b',`e'[1,4],`e'[1,1],`e'[1,2],`e'[1,3],`=`e'[1,3]+1' //elasticity
			loc names `names' alphaH barK mu epsilon epsilon_plus1
			loc coleq `coleq' tcr tcr tcr tcr
		}
		if "`kink'"!="" {
			mat `b'=`b',`=ln(`=eresp'/`cutoff'+1)/(ln(1-`t0_kink')-ln(1-`t1_kink'))'
			loc names `names' elasticity
			loc coleq `coleq' kink
		}	
		if "`notch'"!="" {
			mata: e=notch(`t0_notch',`t1_notch',`deltaT_notch',`cutoff',`=eresp',`init')
			mata: st_matrix("`e'",e)
			if `e'[1,2]!=0 {
				noi di as error "Error code `errorcode' during numeric optimization using notch(). See help mata optimize."
				noi di as error "Elasticity estimates not reported"
				mat `b'=`b',. //elasticity
			}
			else mat `b'=`b',`e'[1,1] //elasticity
			loc names `names' elasticity
			loc coleq `coleq' notch
		}


		mat colnames `b'=`names'
		mat coleq `b'=`coleq'

		eret post `b', esample(`touse') obs(`N')
		ereturn matrix polynomial=`polynomials'
		ereturn scalar bandwidth=`bw'
		ereturn scalar cutoff=`cutoff'
		ereturn scalar lower_limit=`zL'
		ereturn scalar upper_limit=`zH'
		if "`localbw'"=="rdbw" {
			mat rownames `localbws'=`yvars' `characterize'
			ereturn matrix localbws = `localbws'
		}
		//ereturn scalar min=`lo'
		//ereturn scalar max=`hi'
		
		ereturn local binname `varlist'
		ereturn local indepvars `xvars'
		mat colnames `table'=`colfreq'
		ereturn matrix table=`table'
		if "`adjust'"!="none"&"`adjust'"!="" {
			mat colnames `adj_table'=`adjnames'
			ereturn matrix adj_table=`adj_table'
			ereturn local adjustment="`adjust'"
			}
		if "`xvars'"!="" ereturn matrix xtypes=`xtypes'
		ereturn local cmdname "rfbunch"
		noi eret di

	}
end

///MATA FUNCTIONS
mata
mata clear  

void evaltcr(todo,p, t , tau , eta , r , eresp,shift, v, g, H)
	 {
	 	
		Ktau=p[1]
		mu=p[2]
		e=p[3]
		alpha=p[4]
		
			
v=		((1-t)*alpha^(1/(e+1))*Ktau^(-1/(e+1))-r+(r/(mu+1))*Ktau^(-mu/(mu+1))*(tau*((mu+1)/mu))^(mu/(mu+1)))^2 \ //FOC for capital, marginal buncher
		((1-(1-eta)*t)^(e+1)*(r*(1-((r^mu)/((mu+1)))))^(-e)*alpha/e-(1-t)*(((e+1)/e)*alpha^(1/(e+1))*Ktau^(e/(e+1))-tau)+r*Ktau-r*(tau*((mu+1)/mu))^(mu/(mu+1))*Ktau^(1/(mu+1)))^2 \ //indifference condition, marginal buncher
		/*(log(mu+1)-log(mu)+log(tau+eresp)+(mu-e)*log(1-t)+(e-mu)*log(r)+(e+1)*log(1-((r^mu)/((mu+1)))*(1-t)^(-mu))-log(alpha))^2 \ //productivity of marginal buncher
		(-log(shift)+(e-mu)*log(1-t)+(e+1)*(-log(1-(1-eta)*t)-log((1-((r^mu)/(mu+1))*(1-t)^(-mu)))+log(1-r^mu/(mu+1))))^2 //ratio of optimal debt
		*/ (((mu+1)/mu)*(tau+eresp)*(1-t)^(mu-e)*r^(-(mu+1))*(r*(1-((r^mu)/((mu+1)))*(1-t)^(-mu)))^(e+1)-alpha)^2 \
		((((1-((r^mu)/((mu+1)))))/(1-((r^mu)/((mu+1)))*(1-t)^(-mu)))^(e+1)*(1-(1-eta)*t)^(-(e+1))*(1-t)^(e-mu)-shift )^2
	}

	function tcr(real scalar t, real scalar tau, real scalar eta, real scalar r, real scalar eresp, real scalar shift, real matrix init)
		{
		mu=init[1]
		e=init[2]
		Ktau=(1-t)^(mu+1)*((mu+1)/mu)*r^(-(mu+1))*tau
		alpha=((mu+1)/mu)*(tau+eresp)*(1-t)^(mu-e)*r^(-(mu+1))*(r*(1-((r^mu)/((mu+1)))*(1-t)^(-mu)))^(e+1)
		S = optimize_init()
		optimize_init_which(S,"min")
		optimize_init_conv_ptol(S, 1e-12)
		optimize_init_conv_vtol(S, 1e-12)
		optimize_init_evaluator(S, &evaltcr())
		optimize_init_evaluatortype(S, "v0")
		optimize_init_argument(S, 1, t)
		optimize_init_argument(S, 2, tau)
		optimize_init_argument(S, 3, eta)
		optimize_init_argument(S, 4, r)
		optimize_init_argument(S, 5, eresp)
		optimize_init_argument(S, 6, shift)
		optimize_init_params(S, (Ktau,mu,e,alpha))
		optimize_init_conv_maxiter(S,10000)
		p=optimize(S)
		return(p,optimize_result_errorcode(S))
		}
	
	
	void evalnotch(todo,e, t0 , t1 , deltaT, cutoff , eresp,v, g, H)
	 {
		v = ( (1/(1+eresp/cutoff))*(1+(deltaT/cutoff)/(1-t0))-(1/(1+1/e))*(1/(1+(eresp/cutoff)))^(1+1/e)-(1/(1+e))*(1-(t1-t0)/(1-t0))^(1+e))^2
	}

	function notch(real scalar t0, real scalar t1, real scalar deltaT, real scalar cutoff, real scalar eresp, real scalar init)
		{
		S = optimize_init()
		optimize_init_which(S,  "min")
		optimize_init_conv_ptol(S, 1e-12)
		optimize_init_conv_vtol(S, 1e-12)
		optimize_init_evaluator(S, &evalnotch())
		optimize_init_evaluatortype(S, "d0")
		optimize_init_argument(S, 1, t0)
		optimize_init_argument(S, 2, t1)
		optimize_init_argument(S, 3, deltaT)
		optimize_init_argument(S, 4, cutoff)
		optimize_init_argument(S, 5, eresp)
		optimize_init_params(S, init)
		optimize_init_conv_maxiter(S,100)
		e=optimize(S)
		return(e,optimize_result_errorcode(S))
		}
		
	function eresp(real scalar B,real scalar tau,real matrix cf, real scalar bw)
		{
		integral=polyinteg(cf,1)
		integral[1]=-polyeval(integral,tau)-B*bw
		roots=polyroots(integral)
		realroots=Re(select(roots, Im(roots):==0))
		out=sort(select(realroots,realroots:>tau)',1)'
		if (cols(out)==0) {
			return(.)
		}
		else return(out[1]-tau)
		}

		
	function shifteval(real matrix X, real scalar zL,real scalar zH,real scalar k,real scalar BM, real scalar bw,real scalar type, real scalar precision, real scalar init,real scalar fill,cutoff,hole)
{
	max=max(X)
	shift=init
	for (i=2;i<=precision;i++) {
		v=1
		while (v>0)	{	
			shift=shift-1/10^(i-1)
			data=fill(X,bw,zL,zH,shift,type,fill,cutoff,hole)
	    
			xbin=J(rows(data[.,1]),1,1)
			
			for (p=1; p<=k; p++) xbin=xbin,data[.,1]:^p
			b=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
			
			if (type==2) 		v=(BM*bw-polyeval(polyinteg(b,1),max/(1+shift))+polyeval(polyinteg(b,1),zL))
			else if (type==3) 	v=(BM*bw-polyeval(polyinteg(b,1),max-log(1+shift))+polyeval(polyinteg(b,1),zL))
			else v=(BM*bw-polyeval(polyinteg(b,1),max)+polyeval(polyinteg(b,1),zL))
			}
		shift=shift+1/10^(i-1)
	}

return(b,shift)
}

function fill(real matrix X,real scalar bw,real scalar zL, real scalar zH, real scalar shift, real scalar type,fill,cutoff,hole) 
	{
		min=min(X)
		max=max(X)
		if (hole==1) zH=min(select(X,X:>cutoff))
		if (type<2&hole==0) {
			bin=ceil((X:-cutoff:-2^-23)/bw):*bw:+cutoff:-bw/2
			y=mm_freq(bin):/(1+(type==1)*shift)
		}	
		else {
			if (type==3) bin=(X:<=cutoff):*(ceil((X:-cutoff:-2^-23):/bw)*bw:+cutoff:-bw/2) :+ ((X:>cutoff):*(floor((X:-zH:+2^-23):/bw):*bw:+zH:-log(1+shift):+bw/2))
			else bin=(X:<=cutoff):*(ceil((X:-cutoff:-2^-23):/bw)*bw:+cutoff:-bw/2) :+ (X:>cutoff):*(floor(((X:-zH:+2^-23):/(1+shift)):/bw):*bw:+zH:/(1+shift):+bw/2)
			y=mm_freq(bin)
		}
		bin=uniqrows(bin)
		
		if (fill==1) {
			if (type<2) 		fullbin=		((zL:-(ceil((zL-min)/bw)::1):*bw) \ (zH:+(0::floor((max-zH)/bw)):*bw)):+bw/2
			else if (type==2) 	fullbin=		((zL:-(ceil((zL-min)/bw)::1):*bw) \ (zH:/(1+shift):+(0::floor(((max-zH)/(1+shift))/bw)):*bw)):+bw/2
			else if (type==3) 	fullbin=		((zL:-(ceil((zL-min)/bw)::1):*bw) \ (-log(1+shift):+zH:+(0::floor(((max-zH))/bw)):*bw)):+bw/2
			if (rows(y)==rows(fullbin)) {
				 fully=y
				}
			else {
				fully=J(rows(fullbin),1,.)
				l=1
				for (j=1;j<=rows(fullbin);j++) {
					if (abs(bin[l]-fullbin[j])<bw/10) {
						fully[j]=y[l]
						l=l+1
					}
					else {
						fully[j]=0
					}
				}
			}
		return(fullbin,fully)
		}
		else return(bin,y)
		}

	

end
