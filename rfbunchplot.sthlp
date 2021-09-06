{smcl}
{cmd:help rfbunchplot}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:rfbunchplot} {hline 2}}Plotting of reduced form bunching estimates from rfbunch.
{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 10 15 2}
{cmd:rfbunchplot} [{varlist}] [{cmd:,} graph_opts(string) parameters(string) noci nostar adjust limit(numlist) weight]


{pstd}
{cmd:rfbunchplot} plots bunching plots after rfbunch estimation for the running variable (if nothing or the name of the running variable is specified as {varlist}) or for an alternative response variable if that is specified.


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt parameters(string)}} Add specified parameters to figure.{p_end}
{synopt:{opt noci}} Do not plot confidence intervals.{p_end}
{synopt:{opt nostar}} Do not add stars for statistical significance (default * 0.1 ** 0.05 *** 0.01).{p_end}
{synopt:{opt adjust}} Plot adjusted data, if used. Only relevant for running variable.{p_end}
{synopt:{opt limit(numlist)}} Only plot values of z between the two numbers in limit().{p_end}
{synopt:{opt weight}} Use bin size to show relative size of bins. Only relevant for alternative response variables..{p_end}
{synopt:{opt graph_opts(string)}} any other option suitable for for twoway.{p_end}
{synoptline}

{marker Author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}Statistics Norway{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}

{marker also_see}{...}
{title:Also see}

{p 7 14 2}
{helpb rfbunch} for help on the main command rfbunch.{p_end}
