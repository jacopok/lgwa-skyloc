---
# compile with:
# pandoc -f markdown-implicit_figures --citeproc results.md -o results.pdf

title: "LGWA sky localization: BNS and BWD"
author: Jacopo Tissino
date: July 2022
geometry: "left=3.5cm,right=3.5cm,top=3cm,bottom=3cm"
output: pdf_document
colorlinks: true
linkcolor: blue
link-citations: true
bibliography: [LGWA.bib]
abstract: |
  What are the detection prospects for the Lunar Gravitational Wave Array?
  What are its sky localization capabilities, specifically regarding 
  neutron stars and white dwarfs?


hyperrefoptions:
- linktoc=all
- linkcolor=blue

---


# LGWA

The Lunar Gravitational Wave Antenna [@harmsLunarGravitationalwaveAntenna2021] is 
a proposed Moon-based gravitational wave detector, most sensitive in the deciHertz band.

![BNS SNR evolution](figs/ET_LGWA_BNS.pdf)

## White dwarfs

White dwarfs are electron-degenerate stars.
They are more compact than main sequence stars, but still much less so than neutron stars
or black holes: typical compactesses $C = R / M$ are on the order of $10^{4}$.
Because of this, they merge at much lower frequencies than compact objects, despite 
being less massive than them.

### BWD merger frequency

How does the merger frequency for white dwarf binaries depend on mass?

Let us work in the Newtonian approximation; then, we know that Kepler's third law holds, 
$$ a^3 f^2 = \frac{G (m_1 + m_2 )}{4 \pi^2}
\,,
$$
and therefore we will have a merger rougly when twice the minor axis of the orbit, $b = a \sqrt{1 - e^2}$,
equals the sum of the radii of the objects, which we assume depend only the masses: 
$$ b = a \sqrt{1 - e^2} = \frac{r(m_1 ) + r(m_2 ) }{2}
\,.
$$

Tidal disruption may occur before this moment, but here we are just trying
to get a scaling relation.

This allows us to determine the orbital frequency $f$ associated to the merger: 
$$ f = \sqrt{\frac{G (m_1 + m_2 )}{4 \pi^2}} \left(\frac{r(m_1 ) + r(m_2 ) }{2 \sqrt{1 - e^2}}\right)^{-3/2} \sim m^{1/2} r^{-3/2}
\,.
$$

The GW emission frequency will be twice this one: $f _{\text{GW}} = 2 f$.

Even though this is only a Newtonian expression, we recover the correct scaling for 
black hole mergers: the Schwarzshild radius scales with $r \propto m$, so $f \propto m^{-1}$.

How does this play out for white dwarf binaries? Their radii depend on their masses in a 
much weaker way, which can be approximated by [@maselliBinaryWhiteDwarfs2020]:
$$ r (m) \approx 0.012 R_{\odot} \sqrt{ 
    \left(\frac{m}{1.44M_{\odot}}\right)^{-2/3} -
    \left(\frac{m}{1.44M_{\odot}}\right)^{2/3}
}
$$
which, as expected, approaches zero for $M \to M _{\text{Chandrasekhar}}$.
For typical white dwarf masses, $M \sim 0.6M_{\odot}$, the first term is more relevant than the second, 
so roughly $r \propto m^{-1/3}$.

Not only is this a much weaker dependence than that of black holes, it's in the opposite direction!
It yields $f \propto m$.
With the complete formula, the dependence is even steeper. 

![Merger frequency as a function of WD mass.](figs/bwd_merger_frequencies_by_mass.pdf)

These frequencies lie squarely inside the LGWA band, which means that the
high-frequency part of this band is not useful to detect them.
The "best" BWD mergers to look out for are therefore the high-mass ones, which 
merge at higher frequency.

## Adding `mlgw_bns` to `GWFish`

Do the SNRs match (with `gwfish_TaylorF2`, or with `lalsim_IMRPhenomD_NRTidal`? Yes! 
Do the PE errors match? No! 

Reference check: comparing `gwfish_TaylorF2` with `lalsim_IMRPhenomD`, almost 
everything agrees to within a few%, with the exceptions of the masses going up to 
0.6 in log-error.

Comparing `mlgw_bns` and `lalsim_IMRPhenomD_NRTidal`: 
things depending on amplitude are kind of OK; masses, phase, (also psi a bit)
are way messed up.

Typically, it seems like errors on the masses are _underestimated_ by $e^{1 \div 4}$.

## Bibliography