{smcl}
{cmd:help rfbunch}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:rfbunch} {hline 2}}Reduced form and structural bunching estimation with multiple response variables.
{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 10 15 2}
{cmd:rfbunch} {depvar} [{indepvars}] {cmd:,} cutoff(real) bw(real) [limits(numlist) notch(numlist) kink(numlist) tcr(numlist) pol(numlist) adjust(string) constant nofill]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Reduced form bunching}
{synopt:{opt bw(real)}} Specifies the bandwidth of the bins; required.{p_end}
{synopt:{opt cut:off(real)}} Specifies the cutoff z*; required.{p_end}
{synopt:{opt lim:its(numlist)}} specify L and H, the number of bins excluded below and above the cutoff. Default is 1 and 0.{p_end}
{synopt:{opt pol:ynomial(numlist)}}specify the degree of the polynomial for the density of the running variable (default=7) and other response variables (default=1). {p_end}
{synopt:{opt adj:ust(numlist)}}Iteratively adjusts the running variable above the cutoff to account for distortions. Allows "none" (default), "y" (Chetty et al., 2011) and "x" (Andresen and Thorvaldsen, 2021).{p_end}
{synopt:{opt constant}}uses the constant approximation to the density in the bunching region (Kleven and Waseem, 2013) rather than the estimated polynomial.{p_end}
{synopt:{opt nofill}}does not fill empty bins, unlike the default which is to use frequency of 0 for any bins where there are no agents.{p_end}

{syntab:Structural}
{synopt:{opt kink(numlist)}}Report the structural elasticity estimates from a tax kink specified by two numbers t0, t1. See Saez 2010, Chetty et al. 2011.{p_end}
{synopt:{opt notch(numlist)}}Report the structural elasticity estimates from a tax notch specified by two or three numbers t0, t1, deltaT. See Kleven and Waseem, 2013.{p_end}
{synopt:{opt tcr(numlist)}}Report structural elasticities of the thin capitalization rule of Andresen and Thorvaldsen (2021).{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rfbunch} calculates reduced form bunching estimates for the main running variable {depvar} and alternative response variables {indepvars} using polynomials.

{marker remarks}{...}
{title:Remarks}

{pstd}
rfbunch first estimates polynomial density estimates of the counterfactual density
 of the running variable, excluding the bunching region. If "adjust" is specified,
 rfbunch repeatedly shifts the running variable above the cutoff until the missing
 mass above the cutoff is equal to the excess mass in the bunching region, either 
 by shifting bin counts up ("y", Chetty et al. 2011) or by shifting observations 
 to the right and reconstructing bins ("x", Andresen and Thorvaldsen, 2021).

{pstd}
If alternative endogenous variables are specified, for each of them rfbunch estimates
 a polynomial relationship between this variable and the running variable separately
 below and above the cutoff and then backs out estimates of the mean of this variable 
 a) among nonbuchers in the bunching region, b) among bunchers and c) for bunchers' 
 counterfactual allocation.

{pstd}
rfbunch does not report standard errors, but standard errors can easily be obtained 
by using bootstrap.

{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:rfbunch} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(polynomial)}}degree of polynomial used for main bunching variable.{p_end}
{synopt:{cmd:e(cutoff)}}Value of cutoff threshold z*{p_end}
{synopt:{cmd:e(lower_limit)}}lower limit of the bunching region{p_end}
{synopt:{cmd:e(upper_limit)}}upper limit of the bunching region.{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rfbunch}{p_end}
{synopt:{cmd:e(binname)}}name of the bunching variable{p_end}
{synopt:{cmd:e(indepvars)}}names of other response variables.{p_end}
{synopt:{cmd:e(properties)}}{cmd:b}; {cmd:V}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(table)}}Table of frequency and bin midpoints in all bins in data, as well as means of all other response variables in those bins.{p_end}
{synopt:{cmd:e(adj_freq)}}adjusted frequencies and bin mindpoints if using adjust(){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{pstd}
Note: {cmd:bootstrap} stores its own output in {cmd:e()} as well.  See
{helpb bootstrap##saved_results:bootstrap}.


{marker references}{...}
{title:References}

{phang}
Andresen, M. and Thorvaldsen, L. (2021) Estimating the tax elasticity of firm debt: 
What can we learn from bunching?, working paper

{phang}
Chetty, Raj, John N. Friedman, Tore Olsen, and Luigi Pistaferri. 2011. Adjustment Costs, 
Firm Responses, and Micro vs. Macro Labor Supply Elasticities: Evidence from Danish Tax 
Records. Quarterly Journal of Economics 126(2): 749-804.

{phang}
Kleven, H. 2016. Bunching. Annual Review of Economics Vol. 8:435-464

{phang}
Kleven, H. and Waseem, M. (2013) Using notches to uncover optimization frictions and 
structural elasticities: Theory and evidence from Pakistan.

{phang}
Saez, E. 2010. Do Taxpayers bunch at kink points? American Economic Journal: 
Economic Policy 2 (August 2010): 180–212


{title:Thanks for citing rfbunch as follows}

{pstd}
Andresen, Martin E., (2021). "RFBUNCH: Stata module to estimate reduced form 
bunching and structural parameters with multiple response variables.". This 
version VERSION_DATE.{p_end}

{pstd}
where you can check your version date as follows:{p_end}

{phang2}{cmd:. which rfbunch}{p_end}

	 
{marker Author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}Statistics Norway{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}

{pstd}
Thanks to John Friedman for advice. This program is work in progress and is used at your own risk.

{marker also_see}{...}
{title:Also see}

{p 4 14 2}
Development version: net install rfbunch, from("https://raw.githubusercontent.com/martin-andresen/rfbunch/master"){p_end}

{p 7 14 2}
{helpb rfbunchplot} for help on the separate plotting command rfbunchplot, for use after rfbunch estimation.{p_end}
