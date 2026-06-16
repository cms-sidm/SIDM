# Pseudoscalar dark-matter bound-state production cross section

This study computes the LHC ($pp$, $\sqrt{s}=13$ TeV) production cross section of a dark-matter
bound state ("darkonium") that mixes with a CP-odd **pseudoscalar mediator** in a dark-matter
simplified model, as a function of the pseudoscalar mass $m_A$, the dark-matter mass $m_\mathrm{DM}$,
and the dark-sector couplings ($\alpha_D$, $y_\chi$, $y_q$, $\tan\beta$).

The physics, every assumption, and the parameter dependence are explained in the notebook
[`pseudoscalar_xsec.ipynb`](pseudoscalar_xsec.ipynb), which is committed with its outputs so the plots
are visible without re-running.

## Files

| file | what it is |
|---|---|
| `pseudoscalar_xsec.ipynb` | the study notebook — physics narrative, parameter scans, and plots |
| `pseudoscalar_xsec.py` | cross-section engine: mediator width, A–B mass mixing, gluon-fusion and $q\bar q$-fusion cross sections, K-factors, scan helpers |
| `mstwpdf.py` | pure-NumPy port of the official MSTW2008 PDF grid reader + bicubic-spline interpolator |
| `mstw2008nnlo.00.dat` | MSTW2008 NNLO central PDF grid (the input the interpolator reads) |
| `mathematica_reference.json` | ground-truth values from the original Mathematica notebook, used by the validation plots (see below) |
| `mathematica/` | the original Mathematica notebook `PscalarXsec.nb` and the inputs needed to run it (provenance) |

The `mathematica/` subfolder holds the original source this study was ported from: `PscalarXsec.nb`
and the MSTW `mstwpdf` package (`mstw2008code/`). The PDF grids are not shipped inside `mathematica/`:
to re-run `PscalarXsec.nb` (Mathematica only, not available on the LPC), download the MSTW2008 NNLO
central grid from the MSTW website into `mstw2008code/Grids/`. The Python port reads its own copy of
that one grid, `mstw2008nnlo.00.dat`, from the top level of this directory.

## Origin: conversion from the Mathematica notebook

This is a port of an original Mathematica notebook, `PscalarXsec.nb`. The port was done so the
calculation runs in the SIDM Python environment on the LPC, where Mathematica is not available, and so
it can be read, reviewed, and reused like the rest of the analysis code.

The conversion was line-by-line and preserves the original physics exactly:

- **PDFs.** The original uses the MSTW collaboration's `mstwpdf` Mathematica package with the
  MSTW2008 NNLO central grid. `mstwpdf.py` is a faithful translation of that package — the same
  bicubic spline over the $(\log_{10}x, \log_{10}Q^2)$ grid, the same heavy-quark threshold handling,
  and the same low-$x$ / high-$Q^2$ extrapolation. PDG flavour codes are used (gluon $=0$).
- **Cross section.** `pseudoscalar_xsec.py` reproduces the original's mediator total width, the
  $2\times2$ complex A–B mass-mixing matrix, the pseudoscalar triangle loop form factor, the
  gluon-fusion and $q\bar q$-fusion cross sections, and the interpolated QCD K-factor tables. The
  parton-luminosity integral is done with a cached Gauss–Legendre rule in $y=\ln x$.

The GeV$^{-2}\to$ fb conversion follows the original, which uses $\hbar c \approx
2\times10^{-14}\,\mathrm{cm\,GeV}$ (i.e. $4\times10^{11}$ fb per GeV$^{-2}$, $\approx2.7\%$ above the
exact $(\hbar c)^2$); the absolute normalization carries that small offset identically in both
implementations (divide absolute rates by $1.027$ to recover the exact normalization).

## Validation and the reference data

`mathematica_reference.json` holds values evaluated **in the Wolfram kernel from the original
notebook's own functions** — the MSTW PDFs $xf(x,Q)$ for several flavours, and the total cross section
$\sigma(m_B)$ scanned over mass for two mediator masses. It is generated on a machine with Mathematica
(not on the LPC) and committed here so the comparison can be reproduced anywhere, including on the LPC
where Mathematica is unavailable. The notebook's validation section loads this file and overlays the
Python results on it:

- the **PDF interpolation** agrees with the Wolfram kernel to a relative $\sim10^{-14}$ (machine
  precision);
- the **total cross section** agrees to $\sim10^{-7}$ (the residual is this Gauss–Legendre
  luminosity integral versus Mathematica's adaptive `NIntegrate`, not a difference in the physics).

The two implementations are independent, so the agreement is a genuine cross-check. It validates
*reproduction of the original calculation*, not its absolute precision: the rate uses a single PDF
member (MSTW2008 NNLO central) and fixed K-factors, so it carries no PDF or QCD-scale uncertainty (the
gluon-fusion scale uncertainty alone is $\mathcal{O}(10$–$20\%)$, well above the $\approx2.7\%$
normalization offset noted above).

### Regenerating the reference data

On a machine with Mathematica, evaluate the original `PscalarXsec.nb` definitions and print the PDF
values `xf[0, x, Q, f]` and the total cross section over a mass scan, then collect them into the same
JSON schema (`pdf` and `sigma` blocks with the couplings recorded under `meta`). The notebook only
reads the JSON, so the LPC side needs no Mathematica.

## Running

Use the `sidm_venv` kernel ("SIDM (LCG_107 Py3.11)"). The notebook adds the repository root to
`sys.path`, imports the two local modules, reads `mstw2008nnlo.00.dat` from this directory, and needs
only `numpy`, `scipy`, `matplotlib`, `mplhep`, and `pandas`.

## References

The model and inputs are documented inline in the notebook. Key references: the ATLAS/CMS Dark Matter
Forum pseudoscalar simplified model (arXiv:1507.00966); the dark-sector bound-state formalism of
Tsai–Wang–Zhao (arXiv:1511.07433) and Elor–Liu–Slatyer–Soreq (arXiv:1801.07723); the pseudoscalar
gluon-fusion (arXiv:1606.00837) and bottom-fusion (hep-ph/0304035) K-factors; and MSTW2008 PDFs
(arXiv:0901.0002).
