# Angular differences

Use `delta_phi` to compare directions across the -pi/pi boundary:

```python
from sidm.tools.utilities import delta_phi

signed = delta_phi(3.1, -3.1)  # approximately -0.083 radians
separation = abs(signed)       # approximately 0.083 radians
```

The result is `phi1 - phi2` wrapped to `[-pi, pi)`, in radians. Scalars,
NumPy arrays and jagged Awkward arrays are supported; missing values are
preserved. Calculations use float64, and opposite directions map to -pi.

For outer-track distances, `dR_outer` uses the wrapped angle and returns
infinity for missing matches or either track's paired `(-5, -5)` coordinates.

Full-range LJ and DSA angular histograms cover 0 to pi, including exact pi.
Generator-lepton Δφ histograms retain their 0–1 radian zoom. Do not merge older
outputs with different angular axes or definitions into these histograms.
