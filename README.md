# RabiHofstadter-CurrentRouting-Code

MATLAB code package for the manuscript:

**Diagonal-hopping control of superradiant photon-current routing in a two-dimensional Rabi--Hofstadter lattice**

This repository provides compact MATLAB routines for representative semiclassical and finite-resource iPEPS single-point calculations. The code computes the current-ordering diagnostics used in the manuscript:

[
Q_{\rm tri}^{+},\quad S_{\rm diag},\quad A_{\rm diag},\quad P_{\rm tri}.
]

## Contents

* `run_four_representative_points.m`
  Main script. Runs the four representative parameter points used in the manuscript.

* `RabiH_DA_mainFunc.m`
  Main function that organizes the semiclassical and iPEPS calculation blocks.

* `RabiH_DA_semiclassicalBlock.m`
  Self-contained semiclassical four-site coherent-state calculation.

* `RabiH_DA_iPEPSBlock.m`
  Compact finite-resource iPEPS single-point calculation block.

* `RabiH_DA_diagnosticsBlock.m`
  Post-processing block for current diagnostics.

* `backend_RHDC_compact.zip`
  Compact backend functions required by the iPEPS block.

* `results_data_availability_four_points/`
  Representative output data generated from the four-point run.

## Representative points

The default script computes the following four points:

[
(\theta_\Delta/\pi,\lambda)=(0,+0.60),\quad
(0,-0.60),\quad
(+0.25,+0.60),\quad
(+0.25,-0.60).
]

These correspond respectively to the square-dominant, triangular-enhanced, diagonal-selected, and mixed routing examples discussed in the manuscript.

## How to run

Open MATLAB and run:

```matlab
run_four_representative_points
```

The script writes output files to:

```text
results_data_availability_four_points/
```

## Notes

The semiclassical block is self-contained. The iPEPS block uses the compact backend included in `backend_RHDC_compact.zip`. The default iPEPS setting is intended as a representative finite-resource reproducibility run. Users may modify the bond dimension, CTM environment dimension, photon cutoff, and model parameters directly in the main script or function files.

## License

This code is released under the MIT License.
