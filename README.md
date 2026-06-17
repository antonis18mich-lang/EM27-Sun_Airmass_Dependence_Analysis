# EM27-Sun Airmass Dependence Analysis

This repository contains a Python module for the analysis of EM27/SUN greenhouse gas retrievals, with a focus on investigating solar zenith angle (SZA) and airmass dependencies in XCO₂, XCH₄, and XAIR measurements.

The main analysis routine reads EM27/SUN retrieval output files, applies quality filtering, calculates summary statistics, and generates diagnostic plots as a function of time, SZA, and airmass. The code can also compare retrievals before and after the application of airmass-dependent correction factors.

To quantify these dependencies, the software performs linear and quadratic regressions and can automatically save fit parameters, including linear slopes and quadratic curvature terms, for later comparison between stations and measurement periods.

The corrected datasets used in the analysis are produced using a separate BASIC-based processing code developed within the COCCON framework. This correction procedure aims to reduce day-to-day variability and systematic dependencies related to solar geometry. The correction code itself is not included in this repository but can be provided upon request.

An example workflow is included in the `Example_Data` folder. Large datasets, plots, and result files are intentionally excluded from the repository to keep it lightweight and portable.

### Main Python Dependencies

* NumPy
* Pandas
* Matplotlib
* SciPy
* Jupyter / IPython (for notebook execution)

This repository is intended as a reusable tool for evaluating airmass-dependent effects and correction strategies in EM27/SUN observations.

