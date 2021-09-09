*! rfbunchplot version date 20210907
* Author: Martin Eckhoff Andresen
* This program is part of the rfbunch package.
cap prog drop rfbunchplot
	program rfbunchplot
	
	syntax [name], [graph_opts(string) parameters(string) noci nostar adjust limit(numlist min=2 max=2) weight]
	
	quietly {
		if "`=e(cmd)'"!="rfbunch" {
			noi di in red "Estimates in memory not created by rfbunch"
			exit
			}
			
		if "`weight'"!="" loc weight [aw=frequency]
		preserve
		if "`namelist'"=="" loc namelist `=e(binname)'
		else if "`namelist'"!="`e(binname)'" {
			if "`adjust'"!="" {
				noi di as error "Option adjust can only be used when plotting main bunching plots - do not specify a variable name other than the running variable"
				exit
			}
			mat b=e(b)
			cap mat b=b[.,"`namelist':"]
			if _rc!=0 {
				noi di as error "Variable `namelist' not present as additional endogenous variable in estimates in memory."
				exit
			}
			
		}
		
		cap confirm matrix e(V)
		if _rc!=0 {
			noi di as text "No variance-covariance matrix found. Confidence intervals and significance stars not reported"
			loc ci noci
			loc star nostar
		}
		
		clear
		
		tempvar f0 f1 CI_l0 CI_r0 CI_l1 CI_r1 error bin f
		mat `f'=e(table)
		svmat `f', names(eqcol)
		rename _`=e(binname)' `=e(binname)'
		rename _frequency frequency
		
		set obs `=_N+1'
		replace `e(binname)'=`e(cutoff)' in `=_N'
		set obs `=_N+1'
		replace `e(binname)'=`=`=_b[bunching:marginal_response]'+`=e(cutoff)'' in `=_N'
		sort `e(binname)'
				
		if "`namelist'"!="`=e(binname)'" {
			gen above=0
			predict double `f0', xb eq(`namelist')
			replace above=1
			predict double `f1', xb eq(`namelist')
			loc eq `namelist'
		}
		else {
			loc eq counterfactual_frequency
			predict double `f0', xb eq(`eq')
		
		}
		
		if "`ci'"!="noci" {
			cap replace above=0
			predict double `error', stdp eq(`eq')
			gen double `CI_l0'=`f0'-invnormal(0.975)*`error'
			gen double `CI_r0'=`f0'+invnormal(0.975)*`error'
			cap replace above=1
			if  "`namelist'"!="`=e(binname)'" {
				drop `error'
				predict double `error', stdp eq(`eq')
				gen double `CI_l1'=`f1'-invnormal(0.975)*`error'
				gen double `CI_r1'=`f1'+invnormal(0.975)*`error'
				}
			if "`=e(binname)'"=="`namelist'" loc ciplot (rarea `CI_l0' `CI_r0' `=e(binname)', color(gs8%50))
			else loc ciplot (rarea `CI_l0' `CI_r0' `=e(binname)' if `=e(binname)'<=`=`=_b[bunching:marginal_response]'+`=e(cutoff)'', color(gs8%50)) (rarea `CI_l1' `CI_r1' `=e(binname)' if `=e(binname)'>=`=e(cutoff)', color(gs8%50))		
			}
				
		foreach param in `parameters' {
			if !inlist("`param'","shift","number_bunchers","share_sample","normalized_bunching","excess_mass","marginal_response","mean_bunchers") {
				noi di as error "List only parameters in  the bunching equation in parameters()."
				exit 301
			} 
			loc eq bunching
			if "`star'"!="nostar" {
				test [`eq']:`param'
				if r(p)<0.01 loc star`param' ***
				else if r(p)<0.05 loc star`param' **
				else if r(p)<0.1 loc star`param' *
			}
			local `param': di %3.2f _b[`eq':`param']
			if "`param'"=="shift" loc `param' `""required shift: ``param''`star`param''""'
			else if "`param'"=="number_bunchers" loc `param' `""number of bunchers: ``param''`star`param'' ""'
			else if "`param'"=="share_sample" loc `param' `""share of sample: ``param''`star`param''""'
			else if "`param'"=="normalized_bunching" loc `param' `""normalized bunching: ``param''`star`param''""'
			else if "`param'"=="excess_mass" loc `param' `""excess mass: ``param''`star`param''""'
			else if "`param'"=="marginal_response" loc `param' `""marginal response: ``param''`star`param''""'
		}
		
		if "`adjust'"!="" {
			if "`namelist'"!="`e(binname)'" {
				noi di as error "Option adjust only for use with basic bunch plots - do not specify alternative plotting variable."
				exit
			}
			if "`e(adjustment)'"=="" {
					noi di as error "No adjustment was made to estimates in e(). Do not specify adjust unless estimates used adjustment".
				exit 
			}
			mat `f'=e(adj_freq)
			svmat `f'
			rename `f'1 adj_bin
			rename `f'2 adj_freq
			
			loc adjplot (bar adj_freq adj_bin if adj_bin>`=e(upper_limit)', barwidth(`=e(bandwidth)') color(maroon%50))
			}
			
			if "`limit'"!="" {
				gettoken min max: limit
				drop if !inrange(`=e(binname)',`min',`max')
			}
			
			su `e(binname)'
			loc xmin=r(min)
			loc xmax=r(max)
			if "`adjust'"!="" {
				su adj_bin
				loc xmax=r(max)
			}
				
			if "`namelist'"=="`e(binname)'" {
				loc lines (line `f0' `e(binname)', color(maroon))
				loc background (bar frequency `e(binname)', color(navy%50) barwidth(`=e(bandwidth)')) 
			}
			else {
				loc lines (line `f0' `e(binname)' if `e(binname)'<=`e(lower_limit)', color(maroon)) (line `f0' `e(binname)' if inrange(`e(binname)',`e(lower_limit)',`=`=_b[bunching:marginal_response]'+`=e(cutoff)''), color(maroon) lpattern(dash)) (line `f1' `e(binname)' if `e(binname)'>=`e(cutoff)', color(maroon))
				loc background (scatter `namelist'b `e(binname)' `weight' if !inrange(`e(binname)',`e(lower_limit)',`e(cutoff)'), color(black) msymbol(circle_hollow)) (scatter `namelist'b `e(binname)' `weight' if inrange(`e(binname)',`e(lower_limit)',`e(cutoff)'), color(maroon))
			}
			
			if "`namelist'"!="`e(binname)'" {
				gen x=`=e(cutoff)'-`e(bandwidth)'/4 in 1
				replace x=`=e(cutoff)'+_b[bunching:average_response] in 2
				replace x=_b[bunching:mean_nonbunchers] in 3
				gen y=_b[`namelist'_means:mean_bunchers] in 1
				replace y=_b[`namelist'_means:mean_bunchers_cf] in 2
				replace y=_b[`namelist'_means:mean_nonbunchers] in 3
				if "`ci'"!="noci" {
					mat ci=e(ci_normal)
					foreach lim in ll ul {
						gen `lim' =ci["`lim'","`namelist'_means:mean_bunchers"] in 1
						replace `lim'=ci["`lim'","`namelist'_means:mean_bunchers_cf"] in 2
						replace `lim'=ci["`lim'","`namelist'_means:mean_nonbunchers"] in 3
					}
					loc scatters (rcap ll ul x, color(dkgreen)) (scatter y x, color(dkgreen))	
				}
				else loc scatters (scatter y x, color(dkgreen))	
			
			if "`ci'"!="noci" loc labels label(1 "95% CI") label(3 "mean in bin") label(5 "polynomial fit") label(9 "estimated means") order(3 5 9 1) cols(2)
			else loc labels label(1 "mean in bin") label(3 "polynomial fit") label(6 "estimated means") order(1 3 6) cols(3)
			loc ytitle `namelist'
			}
			else {
				if "`adjust'"=="" {
					if "`ci'"=="noci" loc labels label(1 "observed") label(2 "estimated counterfactual") order(1 2)
					else loc labels label(1 "95% CI") label(2 "observed") label(3 "estimated counterfactual") order(2 3 1) cols(3)
				}
			else {
				if "`ci'"=="noci" loc labels label(1 "observed") label(2 "adjusted") label(3 "estimated counterfactual") order(1 2 3) cols(3)
				else loc labels label(1 "95% CI") label(2 "observed") label(3 "adjusted") label(4 "estimated counterfactual") order(2 3 4 1) cols(2)
			}
				loc ytitle frequency
			}
			
			if "`namelist'"=="`=e(binname)'" {
				if "`ci'"=="noci" su frequency
				else su `CI_r0'
				loc ymax=r(max)
			}
			else {
				if "`ci'"=="noci" su `=e(binname)'
				else su y
				loc ymax=r(max)
			}

			twoway 	`ciplot' ///
				`background' `adjplot' ///
				`lines'  ///
				`scatters' ///
				, xline(`=e(upper_limit)', lpattern(dash) lcolor(black)) xline(`=e(cutoff)', lpattern(dash) lcolor(maroon)) xline(`=e(lower_limit)', lpattern(dash) lcolor(black)) xline(`=`=_b[bunching:marginal_response]'+`=e(cutoff)'', lpattern(dash) lcolor(navy)) xscale(range(`min' `max')) text(`ymax' `xmax' `shift' `number_bunchers' `share_sample' `normalized_bunching' `excess_mass' `marginal_response' `production' `capital', placement(west) justification(left) size(small)) graphregion(fcolor(white) lcolor(white)) plotregion(lcolor(black)) bgcolor(white) ytitle(`ytitle') legend(`labels') `graph_opts'


		restore
	}
	
	end
