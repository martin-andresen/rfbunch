*! rfbunchplot version date 20220111
* Author: Martin Eckhoff Andresen
* This program is part of the rfbunch package.
cap prog drop rfbunchplot
	program rfbunchplot
	
	syntax [name], [graph_opts(string) /*parameters(string)*/ noci nostar adjust limit(numlist min=2 max=2) weight nomeans]
	quietly {
		if "`=e(cmdname)'"!="rfbunch" {
			noi di in red "Estimates in memory not created by rfbunch"
			exit
			}
		
		if "`weight'"!="" loc weight [aw=frequency]
		preserve
		if "`namelist'"=="" loc namelist `=e(binname)'
		else if "`namelist'"!="`e(binname)'" {
			mat b=e(b)
			cap mat b=b[.,"f0_`namelist':"]
			if _rc!=0 {
				noi di as error "Variable `namelist' not present as additional endogenous variable in estimates in memory."
				exit
			}
			tempname xtypes
			mat `xtypes'=e(xtypes)
			loc xtype=`xtypes'["`namelist'",1]	
		}
		
		if "`means'"!="nomeans"&"`namelist'"!="`e(binname)'" {
			cap confirm number `=_b[`namelist'_effects:bunchers_mean]'
			if _rc!=0 {
				noi di as text "Note: No means for bunchers and non-bunchers found in estimates. These are not plotted."
				loc means nomeans
			}
		}
		
		cap confirm matrix e(V)
		if _rc!=0 {
			noi di as text "No variance-covariance matrix found. Confidence intervals and significance stars not reported."
			loc ci noci
			loc star nostar
		}
		
		clear
		
		tempvar h0 f0 f1 h0_l h0_u f0_l f0_u f1_l f1_u freq meanmat x cimat stdp
		mat `freq'= e(table)
		svmat `freq', names(col)
		
		
		
		if "`limit'"!="" {
			gettoken min max: limit
			drop if `e(binname)'<`min'|`e(binname)'>`max'
		}
		su `e(binname)'
		loc xmin=r(min)
		loc xmax=r(max)
		
		loc plus=0
		if "`namelist'"=="`=e(binname)'" {
			mat `h0'=e(b)
			mat `h0'=`h0'[1,"counterfactual_frequency:"]
			mat score `h0'=`h0'
			
			
			if "`ci'"!="noci" {
				predict double `stdp', stdp
				gen `h0_l'=`h0'-invnormal(0.975)*`stdp'
				gen `h0_u'=`h0'+invnormal(0.975)*`stdp'
				loc ciplot (rarea `h0_l' `h0_u' `e(binname)', color(gs6%50))
				loc cilabno=2
				loc ++plus
				loc cilab label(`cilabno' "95% CI")
			}
		}
			
		loc zL=e(lower_limit)
		loc zH=e(upper_limit)
		

		if "`e(binname)'"!="`namelist'" {
			forvalues t=0/1 {
				mat `f`t''=e(b)
				mat `f`t''=`f`t''[1,"f`t'_`namelist':"]
				mat score `f`t''=`f`t''
			}
			if "`means'"!="nomeans" {
				mat `meanmat'=e(b)
				mat `meanmat'=`meanmat'[1,"`namelist'_effects:"]'
				svmat `meanmat'
				loc ytitle `namelist'
				if `xtype'==1|`xtype'==2 loc add=1
				else loc add=0
				gen `x'=_b[bunching:mean_h0L] in `=1+`add''
				replace `x'=_b[bunching:mean_h0H] in `=2+`add''
				replace `x'=`zL'+(`=e(cutoff)'-`zL')*0.9 in `=3+`add''
				replace `x'=`zH'-(`zH'-`=e(cutoff)')*0.9 in `4+`add''
				loc meanscatter (scatter `meanmat'1 `x', color(dkgreen) msymbol(circle))
				loc ++plus
				loc meanlabno=3
				loc meanlab label(`meanlabno' "estimated means")
				if "`ci'"!="noci"{
					mat `cimat'=e(ci_normal)
					mat `cimat'=`cimat'[.,"`namelist'_effects:"]'
					svmat `cimat'
					loc ciplot (rcap `cimat'1 `cimat'2 `x' in 1/4, color(dkgreen))
					loc cilabno=4
					loc cilab label(`cilabno' "95% CI") 
					loc ++plus
					}
				}
			}
		else loc ytitle frequency
		
		if "`namelist'"=="`e(binname)'" {
			loc background (bar freq `e(binname)', barwidth(`=e(bandwidth)') color(navy))
			loc lines (line `h0' `e(binname)' if `e(binname)'<=`zL', color(maroon)) (line `h0' `e(binname)' if `e(binname)'>`zH', color(maroon)) (line `h0' `e(binname)' if `e(binname)'>`zL'&`e(binname)'<=`zH', color(maroon) lpattern(dash))
			loc labels label(1 "observed") label(`=2+`plus'' "estimated counterfactual") `cilab' order(1 `=2+`plus'' `cilabno') cols(`=2+`plus'')
		}
		else {
			loc background (scatter `namelist' `e(binname)' if `e(binname)'<=`zL'|`e(binname)'>`zH', color(black) msymbol(circle_hollow)) (scatter `namelist' `e(binname)' if `e(binname)'<=`zH'&`e(binname)'>`zL', color(maroon) msymbol(circle)) `meanscatter'
			loc lines (line `f0' `e(binname)' if `e(binname)'<=`zL', color(maroon)) (line `f0' `e(binname)' if `e(binname)'>`zL', color(maroon) lpattern(dash))  (line `f1' `e(binname)' if `e(binname)'>`zH', color(navy)) (line `f1' `e(binname)' if `e(binname)'<=`zH', color(navy) lpattern(dash))
			loc labels label(1 "mean in bin") label(`=3+`plus'' "polynomial fit") `meanlab' `cilab' order(1 `=3+`plus'' `meanlabno' `cilabno') cols(`=2+`plus'')
		}
			
		twoway `background' `ciplot' `lines', xline(`zH', lpattern(dash) lcolor(black)) xline(`zL', lpattern(dash) lcolor(black)) xline(`=e(cutoff)', lpattern(dash) lcolor(maroon)) graphregion(color(white)) plotregion(lcolor(black)) ytitle("`ytitle'") legend(`labels') `graph_opts' 	 		

		restore
	}
	
	end
