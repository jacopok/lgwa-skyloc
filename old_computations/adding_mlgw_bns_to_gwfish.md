## Adding `mlgw_bns` to `GWFish`

Do the SNRs match (with `gwfish_TaylorF2`, or with `lalsim_IMRPhenomD_NRTidal`? Yes! 
Do the PE errors match? No! 

Reference check: comparing `gwfish_TaylorF2` with `lalsim_IMRPhenomD`, almost 
everything agrees to within a few%, with the exceptions of the masses going up to 
0.6 in log-error.

Comparing `mlgw_bns` and `lalsim_IMRPhenomD_NRTidal`: 
things depending on amplitude are kind of OK; masses, phase, (also psi a bit)
are way messed up.

This is _not_ an issue with the low-frequency extension, since
even if we exclude it 

Typically, it seems like errors on the masses are _underestimated_ by $e^{1 \div 4}$,
while the error on the phase is way _overestimated_ by .
