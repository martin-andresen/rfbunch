*! rfbunchdata version date 20210914
* Author: Martin Eckhoff Andresen
* This program is part of the rfbunch package.
cap prog drop rfbunchdata
	program rfbunchdata, rclass
	syntax, [eta(real 0) t(real 0.2) cut(real 0.1) e(real 0.2) mu(real 0.2) r(real 0.2) obs(integer 5000) dist(string) seed(numlist min=1 max=1 integer) nodrop nonoise]
		
	if "`seed'"!="" set seed `seed'
	if "`dist'"=="" loc dist rbeta(2,5)
	
	foreach param in eta t cut e mu r {
		if ``param''<0 {
			noi di as error "`param' must be >=0. You specified ``param''"
			exit 301
		}
		if inlist("`param'","`t'","`eta'")&``param''>1 {
			noi di as error "`param' must be <1. You specified ``param''"
			exit 301
			}
		}
	qui {
		loc Q0=`r'*(1-((`r'^`mu')/(`mu'+1))*(1-`t')^(-`mu'))
		loc Q1=`r'*(1-((`r'^`mu')/(`mu'+1)))
		loc alphaL=`cut'*((`mu'+1)/`mu')*`r'^(-(`mu'+1))*(1-`t')^(`mu'-`e')*`Q0'^(`e'+1)
		//noi di `alphaL'
		
		if `eta'>0 loc alphatilde=((`cut'/`eta')*`e'/(`e'+1))*`Q0'^(`e')*(1-`t')^(-`e')
		//noi di `alphatilde'
		
		mata: st_matrix("alphaH",alphaH(`t',`e',`r',`mu',`cut',`eta',(`=`alphaL'*1.5',`=`Q0'^(-(`e'+1))*(1-`t')^(`e'+1)*`alphaL'*1.5')))
		loc alphaH=alphaH[1,1]
		//noi di `alphaH'
		
		if `eta'>0 {
			if `alphatilde'<`alphaL' {
				noi di in red "Warning: Relevant cutoff is the relative cutoff, not the absolute"
				exit
			}
		}
		
		if `alphaL'>`alphaH' {
			noi di in red "Warning: Parameters are such that the marginal buncher does not exist."
			exit 301
		}
		
		clear
		set obs `obs'
		cap gen alpha=`dist'
		if _rc!=0 {
			noi di in red "distribution provided in dist() unsupported. Provide a random distribution supported by Stata, such as rbeta(2,5)."
			exit
		}
		su alpha
		if r(min)<0 {
			noi di in red "distribution provided in dist() generated negative values for alpha. Use a distribution that generates nonrandom numnbers."
			exit
		}
		count if alpha<`alphaL'
		if `=r(N)'==`=_N' {
			noi di in red "Warning: All agents allocate below cutoff"
			exit 301
		}
		
		if `=r(N)'==0 {
			noi di in red "Warning: All agents allocate at or above cutoff"
			exit 301
		}
		
		loc rep=0
		gen Pi1=(alpha/`e')*`Q1'^(-`e')*(1-(1-`eta')*`t')^(`e'+1)
		gen Pi0=(alpha/`e')*`Q0'^(-`e')*(1-`t')^(`e'+1)
		
		
		gen K0=`Q0'^(-(`e'+1))*(1-`t')^(`e'+1)*alpha
		gen K1=`Q1'^(-(`e'+1))*(1-(1-`eta')*`t')^(`e'+1)*alpha
		
		loc rep=0
		sort alpha
		gen KK=.
		count if alpha<`alphaL'
		qui forvalues i=`=r(N)+1'/`=_N' {
			loc ++rep
			if `rep'==1 mata: Ktau=Ktau(`=alpha[`i']',`t',`e',`r',`mu',`cut',`=K0[`i']')
			else mata: Ktau=Ktau(`=alpha[`i']',`t',`e',`r',`mu',`cut',`=KK[`=`i'-1']')
			mata: st_matrix("Ktau",Ktau)
			if (1-`t')*(((`e'+1)/`e')*`=alpha[`i']'^(1/(`e'+1))*KK[`i']^(`e'/(`e'+1))-`cut')-`r'*KK[`i']+`r'*(`cut'*(`mu'+1)/`mu')^(`mu'/(`mu'+1))*KK[`i']^(1/(`mu'+1))<Pi1[`i'] continue, break
			replace KK=Ktau[1,1] in `i'
			}
		
		gen Pistar=(1-`t')*(((`e'+1)/`e')*alpha^(1/(`e'+1))*KK^(`e'/(`e'+1))-`cut')-`r'*KK+`r'*(`cut'*(`mu'+1)/`mu')^(`mu'/(`mu'+1))*KK^(1/(`mu'+1))

		gen buncher=inrange(alpha,`alphaL',`alphaH')
		
		gen debtcost0=(`mu'/(`mu'+1))*(1-`t')^(`e'-`mu')*`Q0'^(-(`e'+1))*`r'^(`mu'+1)*alpha
		gen debtcost1=(`mu'/(`mu'+1))*(1-(1-`eta')*`t')^(`e'+1)*`Q1'^(-(`e'+1))*`r'^(`mu'+1)*alpha
		
		gen debtratio0 = `r'^`mu'*(1-`t')^(-`mu')
		gen debtratio1 = `r'^`mu'
		gen debtratioK=(`cut'*(`mu'+1)/`mu')^(`mu'/(`mu'+1))*KK^(-`mu'/(`mu'+1))
		
		gen equitycost0=(1-debtratio0)*`r'*K0
		gen equitycost1=(1-debtratio1)*`r'*K1
		gen equitycostK=(1-debtratioK)*`r'*KK
		
		gen equity0=(1-debtratio0)*K0
		gen equity1=(1-debtratio1)*K1
		gen equityK=(1-debtratioK)*KK
		
		gen debt0=debtratio0*K0
		gen debt1=debtratio1*K1
		gen debtK=debtratioK*KK
		
		gen y0=((`e'+1)/`e')*`Q0'^(-`e')*(1-`t')^(`e')*alpha
		gen y1=((`e'+1)/`e')*`Q1'^(-`e')*(1-(1-`eta')*`t')^(`e')*alpha
		gen yK=((`e'+1)/`e')*alpha^(1/(`e'+1))*KK^(`e'/(`e'+1))
		
		foreach var in debtratio equitycost K y equity debt {
			gen `var'=`var'K if KK!=.
			replace `var'=`var'0 if alpha<`alphaL'
			replace `var'=`var'1 if alpha>`alphaH'
		}
		
		gen debtcost = `cut'-2^-23 if KK!=.
		replace debtcost=debtcost0 if alpha<`alphaL'
		replace debtcost=debtcost1 if alpha>`alphaH'
			
		gen tcrcost=max(0,`t'*(debtcost1-`eta'*y1))*(debtcost1>max(`eta'*y1,`cut'))
		
		if "`noise'"!="nonoise" {
			replace debtratio=debtratio+rnormal(0,0.25)
			replace K=K+rnormal(0,1)
			replace debtcost=debtcost-runiform(0,`=0.2*`cut'') if KK!=.
		}

		if "`drop'"!="nodrop" drop KK Pi1 Pi0 alpha buncher debtcost0 debtcost1 debtratio1 debtratio0 Pistar equitycost1 equitycost0 equitycostK debtratioK y0 y1 yK K1 K0 tcrcost equity0 equity1 equityK debt0 debt1 debtK
		
		return scalar marginalresp=`=(`mu'/(`mu'+1))*(1-`t')^(`e'-`mu')*`Q0'^(-(`e'+1))*`r'^(`mu'+1)*`alphaH'-`cut''
		return scalar shift=`Q0'^(-`e'-1)*`Q1'^(`e'+1)*(1-`t')^(`e'-`mu')*(1-(1-`eta')*`t')^(-`e'-1)
		loc Ktaumarg=Ktau[1,1]
		foreach param in alphaH alphaL alphatilde eta t cut e mu r Ktaumarg {
			cap return scalar `param'=``param''
			}
		}
end
	

mata
function evalK(todo,K,alpha,t,e,r,mu,tau,v,g,H) {
	v=((1-t)*alpha^(1/(e+1))*K^(-1/(e+1))-r+(r/(mu+1))*K^(-mu/(mu+1))*(tau*((mu+1)/mu))^(mu/(mu+1)))^2
}
function Ktau(alpha,t,e,r,mu,tau,init) {
		S = optimize_init()
		optimize_init_which(S,"min")
		optimize_init_conv_ptol(S, 1e-12)
		optimize_init_conv_vtol(S, 1e-12)
		optimize_init_evaluator(S, &evalK())
		optimize_init_evaluatortype(S, "v0")
		optimize_init_argument(S, 1, alpha)
		optimize_init_argument(S, 2, t)
		optimize_init_argument(S, 3, e)
		optimize_init_argument(S, 4, r)
		optimize_init_argument(S, 5, mu)
		optimize_init_argument(S, 6, tau)
		optimize_init_params(S, (init))
		optimize_init_conv_maxiter(S,100)
		p=optimize(S)
		return(p)
}

function evalalphaH(todo,p,t,e,r,mu,tau,eta,v,g,H) {
	alpha=p[1]
	K=p[2]
	alphaL=tau*((mu+1)/mu)*r^(-(mu+1))*(1-t)^(mu-e)*(r*(1-((r^mu)/(mu+1))*(1-t)^(-mu)))^(e+1)
	v= 	((1-t)*alpha^(1/(e+1))*K^(-1/(e+1))-r+(r/(mu+1))*K^(-mu/(mu+1))*(tau*((mu+1)/mu))^(mu/(mu+1))+(alpha<=alphaL)*10^200)^2 \
		((1-t)*(((e+1)/e)*alpha^(1/(e+1))*K^(e/(e+1))-tau)-r*K+r*(tau*(mu+1)/mu)^(mu/(mu+1))*K^(1/(mu+1))-(1-(1-eta)*t)^(e+1)*(r-r^(mu+1)/(mu+1))^(-e)*alpha/e +(alpha<=alphaL)*10^200)^2 
}

function alphaH(t,e,r,mu,tau,eta,init) {
		S = optimize_init()
		optimize_init_which(S,"min")
		optimize_init_conv_ptol(S, 1e-12)
		optimize_init_conv_vtol(S, 1e-12)
		optimize_init_evaluator(S, &evalalphaH())
		optimize_init_evaluatortype(S, "v0")
		optimize_init_argument(S, 1, t)
		optimize_init_argument(S, 2, e)
		optimize_init_argument(S, 3, r)
		optimize_init_argument(S, 4, mu)
		optimize_init_argument(S, 5, tau)
		optimize_init_argument(S, 6, eta)
		optimize_init_params(S, (init))
		optimize_init_conv_maxiter(S,100)
		p=optimize(S)
		return(p)
}

end

