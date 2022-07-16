---
# compile with:
# pandoc -f markdown-implicit_figures --citeproc results.md -o results.pdf

title: "LGWA sky localization: BNS and BWD"
author: Jacopo Tissino
date: July 2022
geometry: "left=3.5cm,right=3.5cm,top=2.5cm,bottom=2.5cm"
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

## Matched filter SNRs

A quick rundown of how SNRs are computed for evolving and monochromatic signals:
we are assuming we use a linear filter to extract a signal from Gaussian, uncorrelated
zero-mean, stationary noise.

It then turns out that, when using the optimal filter, the SNR will be 
$$ \text{SNR}^2 
= 4 \int_0^{ \infty } \mathrm{d} f \frac{|\widetilde{h}(f)|^2 }{S_n(f)}
= \int_0^{ \infty } \mathrm{d}\log f \left( \frac{h_c}{h_n}\right)^2
\,,
$$
where 
$$ h_c = 2 f |\widetilde{h}(f)|
\qquad \text{and} \qquad
h_n = \sqrt{f S_n(f)}\,.
$$

This assumes that the whole frequency evolution of the signal occurs while we are observing it: 
this is often not the case for very low-frequency inspirals, which may need thousands of years to evolve.

For them, it's more useful to write $|\widetilde{h}(f)| \sim \delta (f-f_0 ) h_0$, where $h_0$ is the time-domain
amplitude of the waveform, and where we use a delta-like function, $\delta(f- f_0 ) \sim T$ for $f\in [f_0, f_0 + 1/T]$
and zero otherwise. 

Then, we get 
$$ \text{SNR}^2
= 4 \int_0^{ \infty } \mathrm{d}f \frac{ (\delta (f-f_0 ) h_0)^2 }{S_n} = 4 \frac{T h_0^2}{S_n(f_0 )}
= \frac{4 T h_0^2}{h_n^2 / f_0} \overset{!}{=} \frac{h_c^2}{h_n^2}
$$
so in order to plot a point at a certain $h_c(f_0)$ we need to set 
$$ h_c = 2\sqrt{T f_0} h_0
$$

# LGWA

The Lunar Gravitational Wave Antenna [@harmsLunarGravitationalwaveAntenna2021] is 
a proposed Moon-based gravitational wave detector, most sensitive in the deciHertz band.

The basic working principle is to have inertial sensing: 
at a very abstract level, suspend a test mass and measure the ground's displacement
compared to it.

## LGWA response and noise

The answer to "how sensitive will LGWA be?" depends on 
how the Moon responds to GWs, how we can get a readout on this,
and what is the noise, both intrinsic and in the readout.

The response is not really well-known: for sure there are
spheroidal normal modes in the mHz band, which will enhance the 
response in certain cases.

The displacement response due to a single mode is in the form 
$$ \xi _n (f) = \frac{1}{2} h(f) L_n T_n(f)
$$
where $T_n$ is a dimensionless transfer function, while $L_n$ is the effective baseline 
for that mode.
$T_n$ is $\sim 1$ for $f \gg f_n$ (the frequency of that mode), 
it has a peak at $f \sim f_n$, and it goes to zero for $f \ll f_n$.

The effective baseline of the various modes is not well-known, but
it can be assumed to be on the order of 0.6 times the Moon's radius, $10^6 \text{m}$.

This gives us an optimistic view of how the GW strain maps to measurable displacement.

We can then combine this with estimates for readout noise sensitivity in order
to get a tentative noise PSD.

## Possible LGWA targets

The decihertz band is a "link" between the space- and ground-based interferometer bands.

Binaries: 

- Early warning for BNS / BBH
- Local LISA verification binaries, BWDs, WD+star, NS+star for long observation
- SMBH/IMBH mergers

## Early warning for CBCs: can we localize them?

First, we show the time-evolution of the characteristic strain of a BNS (with GW170817's parameters)
and a BBH (with GW150914's parameters) as would be seen by LGWA + ET.

![BNS SNR evolution](figs/ET_LGWA_BNS_170817.pdf)
![BBH SNR evolution](figs/ET_LGWA_BBH_150914.pdf)

Since the BNS is $\sim 10$ times less massive but also $\sim 10$ times closer, the position of the 
characteristic strains of both signals is quite similar; 
however, the $t(f)$ map differs significantly between the two. 

I ran a simulation of $10^4$ GW170818-like BNS mergers, all located at 40 Mpc,
with randomized phase, sky position, time of arrival, $\psi$, orientation.

Around 32% of these were above a SNR=9 threshold, and among those we have qualitatively the
following results:

- the SNRs are not far above 9, reaching only as much as roughly 15;
- the time of merger is quite well estimated, with $\sigma _t < 0.1 \text{s}$ always;
- the masses are quite well-measured (too well measured? $10^{-6} M_{\odot}$ errors...);
- the sky localization is quite bad: only a small fraction fall below 100 square degrees.

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

This procedure also works, at least at an order-of-magnitude level, for NSs and BHs:
the following plot is made assuming that all NSs have a radius of 12km while all BHs are nonspinning.

![Merger frequency as a function of mass for WDs, NS, BHs.](figs/merger_freq_plot_wd_to_bh.pdf)

Local BWD population studies were performed for LISA [@korolObservationallyDrivenGalactic2022],
since local binaries which are far from the merger, with periods of 
hours to days, abound in our galaxy, with distances as low as hundreds of parsecs.

The lowest periods observed are $P \sim 10^{-1} \text{hr}$, 
meaning $f_{gw} \sim 5 \times 10^{-3} \text{Hz}$.

## Bibliography