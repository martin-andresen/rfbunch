*! rfbunchplot version date 20211007
* Author: Martin Eckhoff Andresen
* This program is part of the rfbunch package.
cap prog drop rfbunchplot
	program rfbunchplot
	
	syntax [name], [graph_opts(string) parameters(string) noci nostar adjust limit(numlist min=2 max=2) weight nomeans]
	
	quietly {
		if "`=e(cmdname)'"!="rfbunch" {
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
		
		if strpos("`=e(characterize)'","`namelist'")>0 loc charvar=1
		else loc charvar=0
		
		cap confirm matrix e(V)
		if _rc!=0 {
			noi di as text "No variance-covariance matrix found. Confidence intervals and significance stars not reported."
			loc ci noci
			loc star nostar
		}
		
		clear
		
		cap scalar b=_b[bunching:marginal_response]
		if _rc==0 loc marginalresponse=`=_b[bunching:marginal_response]'
		else loc marginalresponse=0
		
		tempvar f0 f1 CI_l0 CI_r0 CI_l1 CI_r1 error bin f
		mat `f'=e(table)
		svmat `f', names(col)
	
		su `e(binname)'
		loc xmin=r(min)
		loc xmax=r(max)
		su `e(binname)' if `e(binname)'>`e(cutoff)'+1e-23
		loc minabove=r(min)
		

		cap confirm matrix e(adj_freq)
		if _rc==0&(("`namelist'"!="`e(binname)'"&`charvar'==1)|("`namelist'"=="`e(binname)'"&"`adjust'"!="")) {
			mat `f'=e(adj_freq)
			svmat `f', names(col)
			su adj_bin
			if r(max)>`xmax' loc xmax=r(max)
			}
		
		sort `e(binname)'
		
		if "`namelist'"=="`e(binname)'" loc eq counterfactual_frequency
		else loc eq `namelist'
		
		rename `e(binname)' bin

		cap set obs 200
		gen `e(binname)'=`xmin'+_n*(`xmax'-`xmin')/200 in 1/200

		gen above=0
		predict double `f0', eq(`eq')
		if `charvar'==0&"`namelist'"!="`e(binname)'" {
			replace above=1
			predict double `f1', eq(`eq')
		}
		
		if "`ci'"!="noci"  {
			replace above=0			
			predict double `error', stdp eq(`eq')
			gen double `CI_l0'=`f0'-invnormal(0.975)*`error'
			gen double `CI_r0'=`f0'+invnormal(0.975)*`error'
			if `charvar'==0&"`namelist'"!="`e(binname)'" {
				drop `error'
				predict double `error', stdp eq(`eq')
				gen double `CI_l1'=`f1'-invnormal(0.975)*`error'
				gen double `CI_r1'=`f1'+invnormal(0.975)*`error'
				drop `error'
			}
		if "`=e(binname)'"=="`namelist'"|`charvar'==1 loc ciplot (rarea `CI_l0' `CI_r0' `=e(binname)' if `e(binname)'<., color(gs8%50))
		else loc ciplot (rarea `CI_l0' `CI_r0' `=e(binname)' if `=e(binname)'<=`=`marginalresponse'+`=e(cutoff)''&`=e(binname)'<., color(gs8%50)) (rarea `CI_l1' `CI_r1' `=e(binname)' if `=e(binname)'<.&inrange(`=e(binname)',`=e(cutoff)',`xmax'), color(gs8%50))	
		}
		
		replace above=`e(binname)'>`e(cutoff)'
		
		//print parameters					
		foreach param in `parameters' {	
			if "`namelist'"=="`=e(binname)'"&!inlist("`param'","shift","number_bunchers","share_sample","normalized_bunching","excess_mass","marginal_response","average_response","total_response","mean_nonbunchers") {
				noi di as error "List only parameters in  the bunching equation in parameters() when plotting standard bunching plot."
				exit 301
			} 
			else if "`namelist'"!="`=e(binname)'"&!inlist("`param'","mean_bunchers","bunchers_diff","excess_value","mean_nonbunchers","mean_bunchers_cf") {
				noi di as error "List only parameters in  `namelist'_means equation in parameters() when plotting alternative endogenous variables."
				exit 301
				}
			if "`namelist'"=="`=e(binname)'" loc eq bunching
			else loc eq `namelist'
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
			else if "`param'"=="average_response" loc `param' `""average response: ``param''`star`param''""'
			else if "`param'"=="total_response" loc `param' `""total response: ``param''`star`param''""'
			else if "`param'"=="mean_nonbunchers" loc `param' `""mean nonbunchers: ``param''`star`param''""'
			else if "`param'"=="bunchers_diff" loc `param' `""bunchers difference: ``param''`star`param''""'
			else if "`param'"=="mean_bunchers" loc `param' `""mean among bunchers: ``param''`star`param''""'
			else if "`param'"=="excess value" loc `param' `""excess value: ``param''`star`param''""'
			else if "`param'"=="mean_nonbunchers_cf" loc `param' `""mean nonbunchers, cf: ``param''`star`param''""'
		}
		
		//adjust option
		if "`adjust'"!="" {
			if "`namelist'"!="`e(binname)'" {
				noi di as error "Option adjust only for use with basic bunch plots - do not specify alternative plotting variable."
				exit
			}
			if "`e(adjustment)'"=="" {
				noi di as error "No adjustment was made to estimates in e(). Do not specify adjust unless estimates used adjustment".
				exit 
			}
			
			loc adjplot (bar adj_freq adj_bin if adj_bin>`=e(upper_limit)*_b[bunching:shift]', barwidth(`=e(bandwidth)') color(maroon%50))
			}
			
		//limit option
		if "`limit'"!="" {
				gettoken xmin xmax: limit
				replace `e(binname)'=. if !inrange(`=e(binname)',`xmin',`xmax')
				replace bin=. if !inrange(bin,`xmin',`xmax')
				cap replace adj_bin=. if !inrange(adj_bin,`xmin',`xmax')
			}

		//Regular plot
		if "`namelist'"=="`e(binname)'" {

			loc lines (line `f0' `e(binname)' if `e(binname)'<., color(maroon))
			loc background (bar frequency bin if bin<. , color(navy%50) barwidth(`=e(bandwidth)')) 
			loc ytitle frequency
			
			if "`adjust'"=="" {
				if "`ci'"=="noci" loc labels label(1 "observed") label(2 "estimated counterfactual") order(1 2)
				else loc labels label(1 "95% CI") label(2 "observed") label(3 "estimated counterfactual") order(2 3 1) cols(3)
			}
			else {
				if "`ci'"=="noci" loc labels label(1 "observed") label(2 "adjusted") label(3 "estimated counterfactual") order(1 2 3) cols(3)
				else loc labels label(1 "95% CI") label(2 "observed") label(3 "adjusted") label(4 "estimated counterfactual") order(2 3 4 1) cols(2)
			}
			
		}
		
		//alternative outcome/characterize plot
			else  {
				if `charvar'==0 {
					loc lines (line `f0' `e(binname)' if `e(binname)'<`e(lower_limit)', color(maroon)) (line `f0' `e(binname)' if inrange(`e(binname)',`e(lower_limit)',`=`e(cutoff)'+`marginalresponse''), color(maroon) lpattern(dash)) (line `f1' `e(binname)' if `e(binname)'>=`minabove'&`e(binname)'<., color(navy)) (line `f1' `e(binname)' if inrange(`e(binname)',`e(cutoff)',`minabove'), color(navy) lpattern(dash)) 
					loc background (scatter `namelist' bin `weight' if !inrange(bin,`=e(lower_limit)',`e(upper_limit)')&bin<., color(black) msymbol(circle_hollow)) (scatter `namelist' bin `weight' if inrange(bin,`e(lower_limit)',`e(cutoff)'), color(maroon))
				}
				else {
					cap su adj_bin if adj_bin>`e(cutoff)'
					cap loc minabove=r(min)
					loc lines (line `f0' `e(binname)' if `e(binname)'<`e(lower_limit)'&`e(binname)'<., color(maroon)) (line `f0' `e(binname)' if `e(binname)'>=`minabove'&`e(binname)'<., color(maroon))  (line `f0' `e(binname)' if inrange(`e(binname)',`e(lower_limit)',`minabove')&`e(binname)'<., color(maroon) lpattern(dash))
					loc background (scatter `namelist' adj_bin `weight' if !inrange(adj_bin,`e(lower_limit)',`e(upper_limit)')&adj_bin<., color(black) msymbol(circle_hollow)) (scatter `namelist' adj_bin `weight' if inrange(adj_bin,`e(lower_limit)',`e(cutoff)')&`e(binname)'<., color(maroon))
				}
				
				if "`means'"!="nomeans" {
					gen x=`=e(cutoff)'-`e(bandwidth)'/4 in 1
					gen y=_b[`namelist'_means:mean_bunchers] in 1
					cap scalar b=_b[bunching:average_response]
					if _rc==0 {
						replace x=`=e(cutoff)'+_b[bunching:average_response] in 2
						replace y=_b[`namelist'_means:mean_bunchers_cf] in 2
						}
					replace x=_b[bunching:mean_nonbunchers] in 3
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
				}
			
			if "`ci'"!="noci"{
				if `charvar'==1 loc labels label(1 "95% CI") label(3 "mean in bin") label(5 "polynomial fit") label(8 "estimated means") order(3 5 8 1) cols(2)
				else loc labels label(1 "95% CI") label(3 "mean in bin") label(5 "polynomial fit") label(10 "estimated means") order(3 5 10 1) cols(2)
			} 
			else {
				if `charvar'==1 loc labels label(1 "mean in bin") label(3 "polynomial fit") label(6 "estimated means") order(1 3 6) cols(3)
				else loc labels label(1 "mean in bin") label(3 "polynomial fit") label(7 "estimated means") order(1 3 7) cols(3)
			}
			loc ytitle `namelist'
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
			
			if `marginalresponse'>0 loc xlinemarg xline(`=`marginalresponse'+`=e(cutoff)'', lpattern(dash) lcolor(navy))
			twoway 	`ciplot' ///
				`background' `adjplot' ///
				`lines'  ///
				`scatters' ///
				, xline(`=e(upper_limit)', lpattern(dash) lcolor(black)) xline(`=e(cutoff)', lpattern(dash) lcolor(maroon)) xline(`=e(lower_limit)', lpattern(dash) lcolor(black)) `xlinemarg' xscale(range(`min' `max')) text(`ymax' `xmax' `shift' `number_bunchers' `share_sample' `normalized_bunching' `excess_mass' `marginal_response' `average_response' `total_response' `mean_nonbunchers' `production' `capital' `excess_value' `mean_bunchers' `mean_nonbunchers_cf' `bunchers_diff', placement(swest) justification(left) size(small)) graphregion(fcolor(white) lcolor(white)) plotregion(lcolor(black)) bgcolor(white) ytitle(`ytitle') legend(`labels') xtitle("`=e(binname)'") `graph_opts'


		restore
	}
	
	end
