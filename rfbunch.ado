	*! rfbunch version date 20220107
	* Author: Martin Eckhoff Andresen
	* This program is part of the rfbunch package.
	cap prog drop rfbunch
	program rfbunch, eclass sortpreserve
		syntax varlist(min=1) [if] [in],  CUToff(real) bw(real) [ ///
		LIMits(string) ///
		notch(numlist min=2 max=3 >=0) ///
		kink(numlist min=2 max=2 >0) ///
		init(numlist min=1 max=2 >0) ///
		POLynomial(numlist min=0 integer >=0) ///
		localbw(string) ///
		localkernel(string) ///
		ADJust(string) ///
		fill ///
		constant ///
		placebo ///
		xtype(numlist <=3 >=0 integer) ///
		precision(integer 3) ///
		]
		
		quietly {

			//check options
			gettoken varlist xvars: varlist
			
			tempname polynomials
			loc numxvars: word count `xvars'
			loc numpoly: word count `polynomial'
			
			if `numpoly'>`numxvars'+1 {
				noi di as error "Specify no more numbers than the number of variables in xvars()+1."
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
						noi di as error "Specify no more numbers in xtype() than the number of x variables in numlist xtypes()."
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

			//limits option
			if "`limits'"!="" {
				loc numvalues: word count `limits'
				if `numvalues'>2 {
					noi di in red "Limits() can only contain one or two entries, either numbers or iterate".
					exit 301
				}
				gettoken L H: limits
				capture confirm integer number `L'
				if _rc!=0 {
					noi di in red "Limits() can only contain one or two entries, either numbers or iterate".
					exit 301
				}
				capture confirm integer number `H'
				if _rc!=0&"`H'"!=" iterate"&"`H'"!="" {
					noi di in red "Limits() can only contain one or two entries, either numbers or iterate".
					exit 301
				}
				if "`H'"=="" loc H=0
				if "`H'"==" iterate"&"`adjust'"!="" {
					noi di in red "Option adjust() cannot be combined with iterate in limits."
					exit 301
				}
				
			}
			else {
				loc L=1
				loc H=0
			}
			
			
			loc zL=`cutoff'-`L'*`bw'
			if "`H'"!=" iterate" loc zH=`cutoff'+`H'*`bw'
			
			cap which moremata.hlp
			if _rc!=0 {
				noi di in red "moremata needed. Install using "ssc install moremata"".
				exit 301
				}
			
			if "`adjust'"!="" {
				if !inlist("`adjust'","y","x","logx") {
					noi di in red "Option adjust can only take values y (Chetty et al.), x and logx (Andresen and Thorvaldsen)"
					exit 301
					}
				}
			if "`adjust'"=="x" loc type=2
			else if "`adjust'"=="logx" loc type=3
			else if "`adjust'"=="y" loc type=1
			else loc type=0
			
			if "`fill'"=="fill" loc fill=1
			else loc fill=0
			
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
			
			//check if there are people in either side
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
			count if `varlist'>`zL'&`varlist'<=`cutoff'
			loc Bunchmass=r(N)
			
			loc colfreq `varlist' frequency
			
			su `varlist' if `varlist'>`cutoff'
			if ((r(min)>`cutoff'+`bw')) {
				noi di as text "Note: Hole detected above cutoff. Bins adjusted above cutoff to avoid half-empty bins."
				loc hole=1
				loc minabove=r(min)
			} 
			else loc hole=0
			if "`adjust'"!="x"|`hole'==1|"`H'"==" iterate" {
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
			
			fvexpand `rhsvars'
			loc names `r(varlist)' _cons
			
			if "`zH'"!="" gen `useobs' = `varlist'<=`zL'|`varlist'>`zH'
			
			//Get counterfactual by adjusting values above z* until missing mass=excess mass
			if inlist("`adjust'","x","y","logx") {
				mata: shift=shifteval(st_data(selectindex(st_data(.,"`useobs'")),"`varlist'"),`zL',`zH',`=`polynomials'[1,1]',`BM',`bw',`type',`precision',0,`fill',`cutoff',`hole')
				mata: st_matrix("`b'",shift)
				mata: st_matrix("`adj_table'",fill(st_data(.,"`varlist'"),`bw',`cutoff',`minabove',`=`b'[1,`=colsof(`b')']',`type',0,`cutoff',`hole'))
				mata: h0=shift[1..cols(shift)-1]
				mat `b'=`b'[1,2..`=colsof(`b')-1'],`b'[1,1],`b'[1,`=colsof(`b')']
				local shift=`b'[1,`=colsof(`b')']
				loc names `names' shift
				loc coleq `coleq' bunching
				loc adjnames adj_bin adj_freq
				if `shift'>0 {
					noi di as text "Note: Estimated shift parameter is positive. You might want to think about whether it makes "
					noi di as text "sense that agents exposed to the regime above the threshold should decrease z relative to their"
					noi di as text "counterfactual under the regime below the cutoff."
				}
				}
			
			//get counterfactual by iteratively adjusting upper bound zH untill missing mass==excess mass
			else if "`H'"==" iterate" {
				mata: zH=iterate(st_data(.,"`varlist'"),`bw',`fill',`zL',`minabove',`cutoff',1,`=`polynomials'[1,1]',`precision',`BM')
				mata: st_matrix("`b'",zH)
				mata: h0=zH[1..cols(zH)-1]
				mat `b'=`b'[1,2..`=colsof(`b')-1'],`b'[1,1],`b'[1,`=colsof(`b')']
				local zH=`b'[1,`=colsof(`b')']
				if r(N)<=1 {
					noi di in red "Iterative procedure to find upper bound did not converge before reaching the maximum of z in data."
					exit 301
				}
				loc names `names' upper_bound
				loc coleq `coleq' bunching
				loc hole=1
				loc zH=`cutoff'+ceil((`zH'-`cutoff')/`bw')*`bw'
				su `varlist' if `varlist'>`zH'
				loc minabove=r(min)
			}
			
			//straightforward with no adjustment
			else {
				mata: data=fill(st_data(selectindex(st_data(.,"`useobs'")),"`varlist'"),`bw',`zL',`zH',0,0,`fill',`cutoff',`hole')
				mata: xbin=J(rows(data[.,1]),1,1)
				mata: for (p=1; p<=`=`polynomials'[1,1]'; p++) xbin=xbin,data[.,1]:^p
				mata:h0=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
				mata: st_matrix("`b'",h0)
				mat `b'=`b'[1,2..`=colsof(`b')'],`b'[1,1]
				local shift=0
			}
			if "`adjust'"=="x" mata: h1=h0:/((1+`shift'):^((0::`=`polynomials'[1,1]')'))
			else if "`adjust'"=="logx" mata: h1=h0[1]+`shift',h0[2..`=`polynomials'[1,1]+1']
			else mata: h1=h0
					
			mata: st_matrix("`table'",fill(st_data(.,"`varlist'"),`bw',`cutoff',`cutoff',0,`type',0,`cutoff',0))
			mata: st_numscalar("maniprangecf",(polyeval(polyinteg(h0,1),`zH')-polyeval(polyinteg(h0,1),`zL'))/`bw')
			mata: st_numscalar("excesscf",(polyeval(polyinteg(h0,1),`cutoff')-polyeval(polyinteg(h0,1),`zL'))/`bw')
			mata: st_numscalar("misscf",(polyeval(polyinteg(h1,1),`zH')-polyeval(polyinteg(h1,1),`cutoff'))/`bw')
			mata: st_numscalar("h0tau",(polyeval(h0,`cutoff')))

			loc B=`Bunchmass'-excesscf
			loc fs=`B'/maniprangecf
			loc fs0=`B'/(excesscf+`B')
			loc fs1=`B'/misscf
			mat `b'=`b',`B',`=`B'/`N'',`fs',`fs0',`fs1' //Number of bunchers, bunchers share of sample, normalized bunching, B relative to counterfanctual in missing region, B relative to counterfactual in excess regione
			loc coleq `coleq' bunching bunching bunching bunching bunching
			loc names `names' number_bunchers share_sample share_manipulationregion share_excessregion share_missregion 
			
			if "`constant'"!="constant" {
				mata: meannonbunch=(polyeval(polyinteg((0,h0),1),`cutoff') -polyeval(polyinteg((0,h0),1),`zL'))/(polyeval(polyinteg(h0,1),`cutoff') -polyeval(polyinteg(h0,1),`zL'))
				mata: st_numscalar("meannonbunch",meannonbunch)
				}
			else scalar meannonbunch=(`cutoff'-`zL')/2
			
			if `B'<0|"`placebo'"!="" {
				if `B'<0&"`placebo'"=="" {
					noi di as text "Negative estimates of B - no bunching in the bunching region. Marginal response,"
					noi di as text "total response and counterfactual mean among bunchers cannot be calculated."
				}
				mat `b'=`b',meannonbunch
				loc coleq `coleq' bunching
				loc names `names' mean_nonbunchers
				}
				
			else {
				if "`constant'"!="constant" {
					mata: eresp=eresp(`B',`cutoff',h0,`bw')
					mata: st_numscalar("eresp",eresp)
					if eresp==. {
						noi di as error 	"No real roots found above the cutoff to the integral of the counterfactual."
						noi di as error 	"Marginal response and various other parameters cannot be estimated. You could"
						noi di as error		"try using the constant approximation to the density in the bunching region by"
						noi di as error		"specifying the option constant."
						exit 301
					}
					
					mata: meanh0L=(polyeval(polyinteg((0,h0),1),`cutoff') -polyeval(polyinteg((0,h0),1),`zL'))/(polyeval(polyinteg(h0,1),`cutoff') -polyeval(polyinteg(h0,1),`zL'))
					mata: meanh0H=(polyeval(polyinteg((0,h0),1),`zH') -polyeval(polyinteg((0,h0),1),`cutoff'))/(polyeval(polyinteg(h0,1),`zH') -polyeval(polyinteg(h0,1),`cutoff'))
					mata: st_numscalar("meanh0L",meanh0L)
					mata: st_numscalar("meanh0H",meanh0H)
					mata: totalresponse=(1/`bw')*(polyeval(polyinteg((0,h0),1),`cutoff'+eresp) -polyeval(polyinteg((0,h0),1),`cutoff'))-(`cutoff'*`B')
					mata: st_numscalar("totalresponse",totalresponse)
					}
				else {
					tempname predcut
					scalar eresp=(`bw'*`B')/h0tau
					scalar meanbunch=eresp*h0tau-`B'*`cutoff'/2
					scalar totalresponse=eresp*h0tau
					}
				
				loc coleq `coleq' bunching bunching bunching bunching bunching
				loc names `names' mean_h0L mean_h0H marginal_response total_response average_response
				mat `b'=`b',meanh0L,meanh0H,eresp,totalresponse,`=totalresponse/`B''
				}
			
			restore
		
			//ESTIMATE REPONSE ALONG OTHER ENDOGENOUS VARS
			if "`xvars'"!="" {
				if `B'<0&"`placebo'"=="" {
					noi di as text "Mean counterfactual for bunchers and difference between bunchers "
					noi di as text "means and this quantity cannot be calculated for alternative because B<0."
				}
				
				preserve
				drop if !`touse'
				keep `varlist' `xvars' `characteruze'
				gen `useobs'=`varlist'<=`zL'|`varlist'>`zH'
				
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
										
				tempname means integerbin adjustbin adjustz pred_excess pred_missing pred_cf
				tempvar predy f0 f1 f above fabove mean_b_cf
				
				if `type'<2&`hole'==0 gen `bin'=ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2
				else gen `bin'=(`varlist'<=`cutoff')*(ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2) + (`varlist'>`cutoff')*(floor((`varlist'-`cutoff'+2^-23)/`bw')*`bw'+`cutoff'+`bw'/2)			
				sort `bin'
				
				egen `integerbin'=group(`bin')
				
				if `type'==2 {
					replace `bin'=(`varlist'<=`cutoff')*(ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2) + ((`varlist'>`cutoff')*(floor(((`varlist'-`minabove'+2^-23)/(1+`shift'))/`bw')*`bw'+`minabove'/(1+`shift')+`bw'/2))
					gen `adjustz'=`varlist'/(1+`shift'*(`varlist'>`cutoff'))
					}
				else if `type'==3 {
					replace `bin'=(`varlist'<=`cutoff')*(ceil((`varlist'-`cutoff'-2^-23)/`bw')*`bw'+`cutoff'-`bw'/2) + ((`varlist'>`cutoff')*(floor((`varlist'-`minabove'+2^-23)/`bw')*`bw'+`minabove'-log(1+`shift')+`bw'/2))
					gen `adjustz'=`varlist'-log(1+`shift')*(`varlist'>`cutoff')				
					}
				egen `adjustbin'=group(`bin')
				
				gen `above'=`varlist'>`cutoff'
				loc i=0
				foreach var in `xvars' `char' {
					loc ++i
					loc fnames0
					loc fnames1
					
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
						loc znames
						forvalues k=1/`=`polynomials'[`=`i'+1',1]' {
							if inlist("`adjust'","","y")/*|`xtypes'[`i',1]==0*/ loc xvar `varlist'
							else loc xvar `adjustz'
							if `k'==1 {
								loc rhsvars c.`xvar'
								loc znames c.`varlist'
							}
							else {
								loc rhsvars `rhsvars'##c.`xvar'
								loc znames `znames'##c.`varlist'
							}
							if `xtypes'[`i',1]==1 {
									loc rhsvarsnl `rhsvarsnl' {beta`k'}*`xvar'^`k'+
							}
						}
					}
					
					if `xtypes'[`i',1]==0  {
						reg `var' `rhsvars' `localweights'  if `varlist'<=`zL' //no link
						mat `f0'=e(b)
						reg `var' `rhsvars' `localweights'  if `varlist'>`zH'
						mat `f1'=e(b)
					}
						
					else if `xtypes'[`i',1]==1  {
						reg `var' `rhsvars' if `varlist'<`zL'
						tempname init
						mat `init'=e(b)
						loc initials
						forvalues k=1/`=colsof(`init')-1' {
							loc initials `initials' beta`k' `=`init'[1,`k']'
						}
						nl (`var'= (`rhsvarsnl' {beta0})*(1+{gamma}*`above') ) if `useobs', initial(beta0 `=_b[_cons]' `initials' gamma 0)
						mat `f0'=e(b)
						mat `f0'=`f0'[1,1..`=`polynomials'[`=`i'+1',1]+1']
						mata:st_matrix("`f1'",st_matrix("`f0'"):*(1+`=_b[gamma:]'))
						loc gamma=_b[gamma:]
					}
					
					else if `xtypes'[`i',1]==2 {
						reg `var' `rhsvars' 1.`above' `localweights'  if `useobs' //assume x0/x1=constant if y is in logs, x1-x0 = constant if y is in levels
						mat `f0'=e(b)
						mat `f0'=`f0'[1,1..`=`polynomials'[`=`i'+1',1]'],`f0'[1,`=colsof(`f0')']
						mat `f1'=e(b)
						mat `f1'=`f1'[1,1..`=`polynomials'[`=`i'+1',1]'],`f1'[1,`=colsof(`f1')']+`f1'[1,`=colsof(`f1')-1']
					}
					
					else if `xtypes'[`i',1]==3 {		
						reg `var' `rhsvars'  `localweights' if `useobs' //characterize, assume x0==x1
						mat `f0'=e(b)
						mat `f1'=e(b)
					}
					
					mat `b'=`b',`f0',`f1'
					fvexpand `znames'
					local names `names' `r(varlist)' _cons `r(varlist)' _cons
					mat `f0'=`f0'[1,`=colsof(`f0')'],`f0'[1,1..`=`polynomials'[`=`i'+1',1]']
					mat `f1'=`f1'[1,`=colsof(`f1')'],`f1'[1,1..`=`polynomials'[`=`i'+1',1]']
					forvalues j=1/`=(`polynomials'[`=`i'+1',1]+1)' {
						loc coleq `coleq' f0_`var'
					}
					forvalues j=1/`=(`polynomials'[`=`i'+1',1]+1)' {
						loc coleq `coleq' f1_`var'
					}

					mata: st_matrix("`pred_excess'",(polyeval(polyinteg(polymult(h0,st_matrix("`f0'")),1),`cutoff')-polyeval(polyinteg(polymult(h0,st_matrix("`f0'")),1),`zL')) /	(polyeval(polyinteg(h0,1),`cutoff')-	polyeval(polyinteg(h0,1),`zL'))) 
					mata: st_matrix("`pred_missing'",(polyeval(polyinteg(polymult(h1,st_matrix("`f1'")),1),`zH')-polyeval(polyinteg(polymult(h1,st_matrix("`f1'")),1),`cutoff')) :/(polyeval(polyinteg(h1,1),`zH')-polyeval(polyinteg(h1,1),`cutoff')))
					mata: w1=(polyeval(polyinteg(h1,1),`zH')-polyeval(polyinteg(h1,1),`cutoff'))/(polyeval(polyinteg(h1,1),`zH')-polyeval(polyinteg(h1,1),`zL'))
					mata: w0=(polyeval(polyinteg(h0,1),`cutoff')-polyeval(polyinteg(h0,1),`zL'))/(polyeval(polyinteg(h0,1),`zH')-polyeval(polyinteg(h0,1),`zL'))
					mata: st_matrix("`pred_cf'",w0*st_matrix("`pred_excess'")+w1*st_matrix("`pred_missing'"))
				
					su `var' if `varlist'<=`cutoff'&`varlist'>`zL'
					loc bunchers_mean=(r(mean)-`pred_excess'[1,1])/`fs0'+`pred_excess'[1,1]
					loc prediction_error0=(r(mean)-`pred_excess'[1,1])
					su `var' if `varlist'>`cutoff'&`varlist'<=`zH'
					loc bunchers_cf=(`pred_missing'[1,1]-r(mean)*(1-`fs1'))/`fs1'
					loc prediction_error1=(`pred_missing'[1,1]-r(mean))
					su `var' if `varlist'<=`zH'&`varlist'>`zL'
					loc itt=r(mean)-`pred_cf'[1,1]
					
					if `xtypes'[`i',1]==0 {
						tempname teff
						mata: st_matrix("`teff'",polyeval(st_matrix("`f1'"),`cutoff')*(1+`shift')-polyeval(st_matrix("`f0'"),`cutoff'))
						loc teff=`teff'[1,1]
						loc names `names' RD_effect
						loc coleq `coleq' `var'_effects
						mat `b'=`b',`teff'
					}
					
					if `xtypes'[`i',1]==1 {
						loc names `names' proportional_shift
						loc coleq `coleq' `var'_effects
						mat `b'=`b',`gamma'
					}
					else if `xtypes'[`i',1]==2 {
						loc names `names' constant_shift
						loc coleq `coleq' `var'_effects
						mat `b'=`b',_b[1.`above']
					}
					loc names `names' predicted_mean_excess predicted_mean_missing prediction_error0 prediction_error1 bunchers_mean bunchers_cf bunchers_late itt sorting_diff relative_sorting
					loc coleq `coleq' `var'_effects `var'_effects  `var'_effects `var'_effects `var'_effects `var'_effects `var'_effects `var'_effects `var'_effects `var'_effects 
					
					mat `b'=`b',`pred_excess'[1,1],`pred_missing'[1,1], `prediction_error0', `prediction_error1',`bunchers_mean',`bunchers_cf',`=`bunchers_mean'-`bunchers_cf'',`itt',`=`bunchers_cf'-`pred_missing'[1,1]',`=`bunchers_cf'/`pred_missing'[1,1]'						
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
			if "`adjust'"!="" {
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
		data=fill(X,bw,zL,zH,shift,type,fill,cutoff,hole)
		xbin=J(rows(data[.,1]),1,1)
		
		for (p=1; p<=k; p++) xbin=xbin,data[.,1]:^p
		b=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
		
		if (type==2) 		v=(BM*bw-polyeval(polyinteg(b,1),max/(1+shift))+polyeval(polyinteg(b,1),zL))
		else if (type==3) 	v=(BM*bw-polyeval(polyinteg(b,1),max-log(1+shift))+polyeval(polyinteg(b,1),zL))
		else v=(BM*bw-polyeval(polyinteg(b,1),max)+polyeval(polyinteg(b,1),zL))
		if (v<0) {
			neg=-1
		}
		else {
			neg=1
		}
		//for (i=2;i<=precision;i++) {
		while (v*neg>0)	{	
			data=fill(X,bw,zL,zH,shift,type,fill,cutoff,hole)
			xbin=J(rows(data[.,1]),1,1)
			for (p=1; p<=k; p++) xbin=xbin,data[.,1]:^p
			b=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
			
			if (type==2) 		v=(BM*bw-polyeval(polyinteg(b,1),max/(1+shift))+polyeval(polyinteg(b,1),zL))
			else if (type==3) 	v=(BM*bw-polyeval(polyinteg(b,1),max-log(1+shift))+polyeval(polyinteg(b,1),zL))
			else v=(BM*bw-polyeval(polyinteg(b,1),max)+polyeval(polyinteg(b,1),zL))
			shift=shift-neg/10^precision
			}
		
		//}

	return(b,shift)
	}

	function iterate(real matrix X, real scalar bw, real scalar fill, real scalar zL, real scalar zH, real scalar cutoff, real scalar hole, real scalar k,real scalar precision, real scalar BM) 
	{
			max=max(X)
			//for (i=2;i<=precision;i++) {
				v=1
				while ((v>0)&(zH<max)) {
					data=fill(select(X,(X:<=zL) :| (X:>zH)),bw,zL,zH,0,0,fill,cutoff,hole)
					
					xbin=J(rows(data[.,1]),1,1)
					for (p=1; p<=k; p++) xbin=xbin,data[.,1]:^p
					b=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
					v=(BM*bw-polyeval(polyinteg(b,1),max)+polyeval(polyinteg(b,1),zL))
					zH=zH+cutoff/10^precision
				}
			//}
			return(b,zH)
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
