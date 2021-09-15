*! rfbunchdata version date 20210914
* Author: Martin Eckhoff Andresen
* This program is part of the rfbunch package.
cap prog drop rfbunchdata
	program rfbunchdata, rclass
	syntax [name], eta(numlist max=1 min=1 >=0 <=1)	t(numlist max=1 min=1  >=0 <=1) cut(numlist max=1 min=1) e(numlist max=1 min=1  >=0) mu(numlist max=1 min=1  >=0) r(numlist max=1 min=1 >=0) [obs(integer 5000) dist(string)]
	
	if "`namelist'"=="" loc namelist debtcost
	if "`dist'"=="" loc dist rbeta(2,5)
	
	qui {
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
			
		loc Q0=`r'*(1-((`r'^`mu')/(`mu'+1))*(1-`t')^(-`mu'))
		loc Q1=`r'*(1-((`r'^`mu')/(`mu'+1)))
		loc alphaL=`cut'*((`mu'+1)/`mu')*`r'^(-(`mu'+1))*(1-`t')^(`mu'-`e')*`Q0'^(`e'+1)
		
		gen Ktau=.
		count if alpha<`alphaL'
		if `=r(N)'==`=_N' {
			noi di in red "Warning: All agents allocate below cutoff"
			exit
		}
		
		if `=r(N)'==0 {
			noi di in red "Warning: All agents allocate at or above cutoff"
			exit
		}
		
		
		loc rep=0
		gen Pi1=(alpha/`e')*`Q1'^(-`e')*(1-(1-`eta')*`t')^(`e'+1)
		sort alpha
		qui forvalues i=`=r(N)+1'/`=_N' {
			loc ++rep
			loc alphaH=alpha[`i']
			if `rep'==1 mata: Ktau=Ktau(`alphaH',`t',`e',`r',`mu',`cut',`Q0'^(-(`e'+1))*(1-`t')^(`e'+1)*`alphaH')
			else mata: Ktau=Ktau(`alphaH',`t',`e',`r',`mu',`cut',`=Ktau[`=`i'-1']')
			mata: st_matrix("Ktau",Ktau)
			replace Ktau=Ktau[1,1] in `i'
			if (1-`t')*(((`e'+1)/`e')*`alphaH'^(1/(`e'+1))*Ktau[`i']^(`e'/(`e'+1))-`cut')-`r'*Ktau[`i']+`r'*(`cut'*(`mu'+1)/`mu')^(`mu'/(`mu'+1))*Ktau[`i']^(1/(`mu'+1))<Pi1[`i'] continue, break
			}
		
		
		
		if `alphaL'>`alphaH' {
			noi di in red "Warning: Parameters are such that the marginal buncher does not exist."
			exit
		}
		
		gen `namelist' = `cut'-1e-6 if Ktau!=.
		replace `namelist'=(`mu'/(`mu'+1))*(1-`t')^(`e'-`mu')*`Q0'^(-(`e'+1))*`r'^(`mu'+1)*alpha if alpha<`alphaL'
		replace `namelist'=(`mu'/(`mu'+1))*(1-(1-`eta')*`t')^(`e'+1)*`Q1'^(-(`e'+1))*`r'^(`mu'+1)*alpha if alpha>`alphaH'
		
		//drop Ktau Pi1 alpha
		
		foreach param in alphaH alphaL eta t cut e mu {
			return scalar `param'=``param''
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
end

