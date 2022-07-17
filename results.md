---
# compile with:
# pandoc -f markdown-implicit_figures --citeproc results.md --toc -o results.pdf

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



## LGWA concept

The Lunar Gravitational Wave Antenna [@harmsLunarGravitationalwaveAntenna2021] is 
a proposed Moon-based gravitational wave detector, most sensitive in the deciHertz band.

The basic working principle is based on seismometers on the moon, which do inertial sensing: 
at a very abstract level, suspend a test mass and measure the ground's displacement
compared to it.

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

We can then combine this with estimates for readout/environmental noise sensitivity in order
to get a tentative noise ASD: if our distance-measurement ASD is $x \text{m} / \sqrt{\text{Hz}}$,
then we can compute a strain ASD like $10^{-6} x / \sqrt{\text{Hz}}$.

This means that the characteristic noise strain will be $h_n = \sqrt{f S_n} = 10^{-6} x \sqrt{f / \text{Hz}}$.
An estimate for $x$ is $5 \times 10^{-13}$, which means we'd get $h_n \sim 5 \times 10^{-19}$; 


## Possible LGWA targets

LGWA will be extremely useful for lunar science, but here we'll focus on the astrophysics.

The decihertz band is the "missing link" between the space- and ground-based interferometer bands
[@seddaMissingLinkGravitationalwave2020].

There are non-binary detection possibilities, such as stochastic backgrounds [@barishImpactMidbandGravitational2021]
or continuous waves from pulsars [@harmsLunarGravitationalwaveAntenna2021], 
but we will focus on binaries.

A decihertz observatory could:

- Give early warning for BNS / BBH which merge in the LVC band (guaranteed!).
- Detect local LISA verification binaries, BWDs, WD+star, NS+star with long observations (basically guaranteed!).
- Detect WD mergers in the local universe (and inform us: are BWD progenitors for SNIa?).
- Detect IMBH mergers (and inform us: do these "seed SMBHs" exist?).

What is the time left for a binary to merge, if it only loses energy through GWs?  
$$ \tau \approx 157 \text{s} \times 
\left(\frac{f}{20 \text{Hz}}\right)^{-8/3}
\left(\frac{m_1 + m_2 }{2.8 M_{\odot}}\right)^{-5/3}
\left(\frac{\nu }{1/4}\right)^{-1}
$$
where $\nu = \mu / M = q / (1+q)^2$ is the symmetric mass ratio, which 
decreases for more asymmetric systems, ranging from 0 to 1/4.

### Matched filter SNRs

A quick rundown of how SNRs are computed for evolving and monochromatic signals:
we are assuming we use a linear filter to extract a signal from Gaussian, uncorrelated
zero-mean, stationary noise.

It then turns out that, when using the optimal such filter, the SNR will be 
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

### Early warning for CBCs: can we localize them?

"[I]t can be expected that [LGWA] would significantly improve parameter estima-
tion of some GW signals compared to LISA alone" [@harmsLunarGravitationalwaveAntenna2021].
I haven't really tested this assumption in the IMBH case (also because the LISA PSD is 
missing from GWFish...), but I have looked at BNS and BWD.

#### SNR evolution

First, we show the time-evolution of the characteristic strain of a BNS (with GW170817's parameters)
and a BBH (with GW150914's parameters) as would be seen by LGWA + ET.

![BNS SNR evolution](figs/ET_LGWA_BNS_170817.pdf)
![BBH SNR evolution](figs/ET_LGWA_BBH_150914.pdf)
![IMBH BBH SNR evolution](figs/ET_LGWA_BBH_IMBH.pdf)

Since the BNS is $\sim 10$ times less massive but also $\sim 10$ times closer, the position of the 
characteristic strains of both signals is quite similar; 
however, the $t(f)$ map differs significantly between the two. 

Similarly, we show an IMBH merger of $10^4 M_{\odot}$ BHs, at $z = 1$.

#### LGWA sky localization

I ran a simulation of $10^4$ GW170818-like BNS mergers, all located at 40 Mpc,
with randomized phase, sky position, time of arrival, $\psi$, orientation.

Around 32% of these were above a SNR=9 threshold, and among those we have qualitatively the
following results:

- the SNRs are not far above 9, reaching only as much as roughly 15;
- the time of merger is quite well estimated, with $\sigma _t < 0.1 \text{s}$ always;
- the masses are quite well-measured (too well measured? $10^{-6} M_{\odot}$ errors...);
- the sky localization is quite bad: only a tiny fraction fall below 100 square degrees,
    while most are worse-localized than the whole sphere.

![SNR and skyloc evolution: ET + LGWA](figs/skyloc.pdf)

This is relevant for WDs as well: the waveform model is basically the same
at such low frequencies (non-point-mass effects are negligible), 
and the distances are similar to the typical close SNIa.
Masses will typically be lower, meaning lower amplitudes, but somewhat closer 
distances will compensate.
Plus, in the case of a WD merger we will not get the full LGWA-band observation.

However, what we _can_ do is to give a pre-alert for ground-based observatories:
although the sky localization is bad, the time-of-merger determination is great!
Also, a double-degenerate-progenitor SNIa would occur in the optical quite close to 
the end of GW emission, allowing for us to match the EM signal.

Such close SNIa are quite bright: they'd have visible magnitudes of around $13$
at $\sim 40 \text{Mpc}$ distances
(compare to the maximum $\sim 18$ of GW170817's kilonova),
since their absolute magnitude is $\sim -19.3$.

### White dwarfs

White dwarfs are electron-degenerate stars, resulting from the evolution of 
low-mass Main Sequence stars.
They are more compact than main sequence stars, but still much less so than neutron stars
or black holes: typical compactesses $C = R / M$ are on the order of $10^{4}$.
Because of this, they merge at much lower frequencies than compact objects, despite 
being less massive than them.

White dwarf binaries are potential SNIa progenitors,
but this has had no observational confirmation yet from SN observations
(although there is evidence from supernova remnants [@schaeferAbsenceExcompanionStars2012]),
which is likely due to these objects' intrinsic faintness 
[@rebassa-mansergasWhereAreDoubledegenerate2019].

#### BWD merger frequency

How does the merger frequency for white dwarf binaries depend on mass?

Let us work in the Newtonian approximation; then, we know that Kepler's third law holds, 
$$ a^3 f^2 = \frac{G (m_1 + m_2 )}{4 \pi^2}
\,,
$$
and we can expect to have a merger roughly when twice the minor axis of the orbit, $b = a \sqrt{1 - e^2}$,
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

![Merger frequency as a function of both WD components' masses.](figs/bwd_time_to_merger_half_minute.pdf)

![Merger frequency as a function of mass for WDs, NS, BHs.](figs/merger_freq_plot_wd_to_bh.pdf)

These frequencies lie squarely inside the LGWA band, which means that the
high-frequency part of this band is not useful to detect them.
The "best" BWD mergers to look out for are therefore the high-mass ones, which 
merge at higher frequency.

This procedure also works, at least at an order-of-magnitude level, for NSs and BHs:
the following plot is made assuming that all NSs have a radius of 12km while all BHs are nonspinning.

Local BWD population studies were performed for LISA [@korolObservationallyDrivenGalactic2022],
since local binaries which are far from the merger, with periods of 
hours to days, abound in our galaxy, with distances as low as hundreds of parsecs.

The lowest periods observed are $P \sim 10^{-1} \text{hr}$, 
meaning $f_{gw} \sim 5 \times 10^{-3} \text{Hz}$.

## Bibliography