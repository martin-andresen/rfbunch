*! rfbunch version date 20210907
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
	ADJust(string) ///
	nofill ///
	constant ]
	
	quietly {

		//check options
		gettoken varlist yvars: varlist
		if "`pol'"=="" {
			loc pol=7
			loc i=0
			foreach yvar in `yvars' {
				loc ++i
				loc pol`i'=1
				}
			}	
		else {
			gettoken pol poly: pol
			loc i=0
			loc last=`pol'
			foreach yvar in `yvars' {
				loc ++i
				if "`poly'"!="" {
					gettoken pol`i' poly: poly
					loc last=`pol`i''
					}
				else loc pol`i'=`last'
				}
			}
		cap which moremata.hlp 
		if _rc!=0 {
			noi di in red "moremata needed. Install using "ssc install moremata"".
			exit 301
			}
		
		if "`adjust'"!="" {
			if !inlist("`adjust'","y","x","none") {
				noi di in red "Option adjust can only take values y (Chetty et al.), x (Andresen and Thorvaldsen) or none"
				exit 301
				}
			}
		if "`adjust'"=="x" loc type=2
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
					noi di in red "syntax of tcr is tcr(t eta r), where t and eta are numbers between 0 and 1"
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
		
		
		tempvar resid freq0 freq touse useobs bin
		tempname table cutvals cf b adj_freq obsbins shiftmat
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
		
		forvalues i=1/`pol' {
			if "`rhsvars'"=="" loc rhsvars  c.`varlist'
			else loc rhsvars `rhsvars'##c.`varlist'
			mat `cutvals'=nullmat(`cutvals') \ `cutoff'^`i'
			loc coleq `coleq' counterfactual_frequency
		}
	
		loc coleq `coleq' counterfactual_frequency
		mat `cutvals'=1 \ `cutvals'
		
		mata: st_matrix("`table'",fill(st_data(.,"`varlist'"),`bw',`zH',`zH',1,0,`fill',`cutoff'))
		
		//Get countefactual and adjust, if using
		if inlist("`adjust'","x","y") {
			mata: shift=shifteval(st_data(selectindex(st_data(.,"`useobs'")),"`varlist'"),`zL',`zH',`pol',`BM',`bw',`type',10,1,`fill',`cutoff')
			mata: st_matrix("`b'",shift)
			mata: st_matrix("`adj_freq'",fill(st_data(selectindex(st_data(.,"`useobs'")),"`varlist'"),`bw',`zL',`zH',`=`b'[1,`=colsof(`b')']',`type',`fill',`cutoff'))
			mat `cf'=`b'[1,1..`=colsof(`b')-1']
			mat `b'=`b'[1,2..`=colsof(`b')-1'],`b'[1,1],`b'[1,`=colsof(`b')']
			fvexpand `rhsvars'
			loc names `r(varlist)' _cons  shift
			loc coleq `coleq' bunching
			}
		else {
			mata: data=fill(st_data(selectindex(st_data(.,"`useobs'")),"`varlist'"),`bw',`zL',`zH',1,0,`fill',`cutoff')
			mata: xbin=J(rows(data[.,1]),1,1)
			mata: for (p=1; p<=`pol'; p++) xbin=xbin,data[.,1]:^p
			mata: b=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
			mata: st_matrix("`b'",b)
			mat `cf'=`b'
			mat `b'=`b'[1,2..`=colsof(`b')'],`b'[1,1]
			fvexpand `rhsvars'
			loc names `r(varlist)' _cons
		}
		
		tempname freq0b freq0tau
		mata: st_matrix("`freq0b'",(polyeval(polyinteg(st_matrix("`cf'"),1),`zH')-polyeval(polyinteg(st_matrix("`cf'"),1),`zL'))/`bw')
		mata: st_matrix("`freq0tau'",(polyeval(st_matrix("`cf'"),`cutoff')))
		loc B=`=`Bunchmass'-`freq0b'[1,1]'
		mat `b'=`b',`B',`=`B'/`N'',`=`B'/`freq0b'[1,1]',`=`B'/`freq0tau'[1,1]' //Number of bunchers, bunchers share of sample, normalized bunching, excess mass
		loc coleq `coleq' bunching bunching bunching bunching bunching bunching bunching
		loc names `names' number_bunchers share_sample normalized_bunching excess_mass marginal_response average_response mean_nonbunchers
		
		if "`constant'"!="constant" {
			mata: eresp=eresp(`B',`cutoff',st_matrix("`cf'"),`bw')
			mata: st_numscalar("eresp",eresp)
			mata: meanbunch=(polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff'+eresp) -polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff'))/(`bw'*`B')-`cutoff'
			mata: st_numscalar("meanbunch",meanbunch)
			mata: meannonbunch=(polyeval(polyinteg((0,st_matrix("`cf'")),1),`cutoff') -polyeval(polyinteg((0,st_matrix("`cf'")),1),`zL'))/(polyeval(polyinteg((st_matrix("`cf'")),1),`cutoff') -polyeval(polyinteg((st_matrix("`cf'")),1),`zL'))
			mata: st_numscalar("meannonbunch",meannonbunch)
			}
		else {
			tempname predcut
			mat `predcut'=`cf'*`cutvals'
			scalar eresp=(`bw'*`B')/`predcut'[1,1]
			scalar meanbunch=eresp/2
			scalar meannonbunch=(`cutoff'-`zL')/2
			}
		mat `b'=`b',eresp,meanbunch,meannonbunch //response of marginal buncher
		
		restore
	
		//ESTIMATE REPONSE ALONG OTHER ENDOGENOUS VARS
		if "`yvars'"!="" {
			preserve
			drop if !`touse'
			keep `varlist' `yvars'
			
			tempname means integerbin
			tempvar predy f0 f1 f above fabove mean_b_cf
			gen `bin'=ceil((`varlist'-`cutoff')/`bw')*`bw'+`cutoff'-`bw'/2
			sort `bin'
			egen `integerbin'=group(`bin')
			gen `above'=`varlist'>`cutoff'
			loc i=0
			foreach var in `yvars' {
				tempname cutvalsy
				loc ++i
				if `pol`i''>0 {
					forvalues k=1/`pol`i'' {
						if `k'==1 loc rhsvars c.`varlist'
						else loc rhsvars `rhsvars'##c.`varlist'
						mat `cutvalsy'=nullmat(`cutvalsy') \ `cutoff'^`k'
					}
					mat `cutvalsy'=1 \ `cutvalsy'
				}
				else mat `cutvalsy'=1
				
				if `pol`i''>0 reg `var' `rhsvars' 1.`above' 1.`above'#(`rhsvars')  if !inrange(`varlist',`zL',`zH')
				else reg `var' 1.`above' if !inrange(`varlist',`zL',`zH')
				mat `b'=`b',e(b)
				mat `f'=e(b)
				if `pol`i''>0 mat `f'=_b[_cons],`f'[1,1..`pol`i'']
				else mat `f'=_b[_cons]
				mata:mean_nonbunchers=(polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`cutoff')-polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`zL'))/(polyeval(polyinteg(st_matrix("`cf'"),1),`cutoff')-polyeval(polyinteg(st_matrix("`cf'"),1),`zL'))
				mata: st_numscalar("mean_nonbunchers",mean_nonbunchers)
				local colnames: colnames e(b)
				loc names `names' `=subinstr("`colnames'","1.`above'","above",.)'
				forvalues j=1/`=colsof(e(b))' {
					loc coleq `coleq' `var'
				}
				predict double `predy' if inrange(`varlist',`zL',`cutoff'), residuals
				su `predy'
				loc excess=r(sum)
				drop `predy'
				
				mata: mean_counterfactual=(polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`cutoff'+`=eresp')-polyeval(polyinteg(polymult(st_matrix("`cf'"),st_matrix("`f'")),1),`cutoff'))/(polyeval(polyinteg(st_matrix("`cf'"),1),`cutoff'+`=eresp')-polyeval(polyinteg(st_matrix("`cf'"),1),`cutoff'))
				mata: st_numscalar("mean_b_cf",mean_counterfactual)
				mat `b'=`b',`excess',mean_nonbunchers,`=`excess'/`B'+mean_nonbunchers',mean_b_cf,`=`excess'/`B''
				loc names `names' excess_value mean_nonbunchers mean_bunchers mean_bunchers_cf bunchers_diff
				loc coleq `coleq' `var'_means `var'_means `var'_means `var'_means `var'_means
				
				reg `var' ibn.`integerbin', nocons
				mat `means'=e(b)
				
				mat `table'=`table',`means''
				loc colfreq `colfreq' `var':b
				}
			restore
			}
		
		//ESTIMATE ELASTICITIES
		tempname e
		if "`tcr'"!="" {
			mata: e=tcr(`t_tcr',`cutoff',`eta_tcr',`r_tcr',`=eresp',`=shift',st_matrix("`initmat'"))
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
			mat `b'=`b',`=`=(`B'*`bw')/`freq0tau'[1,1]' / ( `cutoff'*(ln((1-`t0_kink')/(1-`t1_kink'))))'
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
		ereturn scalar polynomial=`pol'
		ereturn scalar bandwidth=`bw'
		ereturn scalar cutoff=`cutoff'
		ereturn scalar lower_limit=`zL'
		ereturn scalar upper_limit=`zH'
		//ereturn scalar min=`lo'
		//ereturn scalar max=`hi'
		ereturn local binname `varlist'
		ereturn local indepvars `yvars'
		ereturn local cmd rfbunch
		mat colnames `table'=`colfreq'
		ereturn matrix table=`table'
		if "`adjust'"!="none"&"`adjust'"!="" {
			ereturn matrix adj_freq=`adj_freq'
			ereturn local adjustment="`adjust'"
			}
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
		return(out[1]-tau)
		}

		
	function shifteval(real matrix X, real scalar zL,real scalar zH,real scalar k,real scalar BM, real scalar bw,real scalar type, real scalar precision, real scalar init,real scalar fill,cutoff)
{
	max=max(X)
	shift=init
	for (i=1;i<=precision;i++) {
		v=1
		while (v>0)	{	
			shift=shift+1/10^(i-1)
			data=fill(X,bw,zL,zH,shift,type,fill,cutoff)
	    
			xbin=J(rows(data[.,1]),1,1)
			
			for (p=1; p<=k; p++) xbin=xbin,data[.,1]:^p
			b=(invsym(quadcross(xbin,xbin))*quadcross(xbin,data[.,2]))'
			
			if (type==2) v=(BM*bw-polyeval(polyinteg(b,1),max*shift)+polyeval(polyinteg(b,1),zL))
			else v=(BM*bw-polyeval(polyinteg(b,1),max)+polyeval(polyinteg(b,1),zL))
			}
		shift=shift-1/10^(i-1)
	}

return(b,shift)
}
function fill(real matrix X,real scalar bw,real scalar zL, real scalar zH, real scalar shift, real scalar type,fill,cutoff) 
	{
		min=min(X)
		max=max(X)
		if (type<2) {
			bin=ceil((X:-cutoff)/bw):*bw:+cutoff:-bw/2
			y=(1+(type==1)*(shift-1)):*mm_freq(bin)
		}
		else {
			bin=(X:<=zH):*(ceil((X:-cutoff)/bw)*bw:+cutoff:-bw/2) :+ ((X:>zH):*(ceil(shift:*(X:-zH):/bw):*bw:+zH*shift:-bw/2))
			y=mm_freq(bin)
		}
		bin=uniqrows(bin)
		if (fill==1) {
			if (type<2) fullbin=	((zL:-(floor((zL-min)/bw)::0):*bw) \ (cutoff:+(1::ceil((max-zH)/bw)):*bw)):-bw/2
			else fullbin=			((zL:-(floor((zL-min)/bw)::0):*bw) \ (shift*cutoff:+(1::ceil(shift*(max-zH)/bw)):*bw)):-bw/2
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
