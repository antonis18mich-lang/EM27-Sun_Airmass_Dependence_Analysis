# EM27-Sun Airmass Dependence Analysis

This repository contains a Python-based analysis framework developed for the investigation of airmass and solar zenith angle (SZA) dependencies in EM27/SUN greenhouse gas retrievals. The primary objective of the project is to evaluate the stability of retrieved column-averaged dry-air mole fractions and to assess the performance of airmass-dependent correction factors applied to the measurements.

The analysis focuses mainly on the atmospheric quantities XCO₂, XCH₄, and XAIR derived from EM27/SUN observations. The software provides tools for reading retrieval output files, filtering measurements, generating diagnostic plots, performing statistical analyses, and extracting parameters that describe the dependence of the retrieved quantities on solar geometry.

The core functionality is implemented in the module:

```
Modules/analyze_em27_dataset.py
```

The module reads combined retrieval output files and computes a variety of diagnostic quantities. Measurements can be filtered according to solar zenith angle and other quality criteria before analysis. The software then calculates summary statistics, including mean values and the number of valid observations, and produces visualizations that help identify systematic dependencies on SZA and airmass.

For each dataset, the code can generate time series plots, scatter plots against solar zenith angle, and scatter plots against airmass. When corrected retrieval products are available, the module allows direct comparison of measurements before and after correction. This makes it possible to evaluate whether the applied correction factors reduce systematic biases and improve the overall stability of the retrieved greenhouse gas columns.

To quantify the observed dependencies, the software performs both linear and quadratic regression analyses. The resulting fit parameters, including the linear slope and quadratic curvature, can be stored automatically in summary tables for later comparison between stations and measurement periods. These outputs are particularly useful for investigating whether a common set of airmass-dependent correction factors is suitable for different locations and instruments.

The framework has been used to analyse measurements from multiple EM27/SUN stations, including Karlsruhe, Izaña, Sodankylä, Rome, and Thessaloniki. The analysis is intended to support studies of instrument performance, retrieval consistency, and the optimization of correction factors used within the COCCON network.

An example workflow is provided in the `Example_Data` directory. The accompanying notebook demonstrates how to import the analysis module, load a retrieval dataset, and perform the complete analysis procedure. The repository therefore serves both as a reusable analysis package and as a reference implementation for future EM27/SUN studies.

Large observational datasets, generated plots, result files, and station-specific archives are intentionally excluded from the repository. Only the analysis code and a minimal example are included in order to keep the repository lightweight and easily portable across different systems.

The code is written in Python and primarily relies on the following scientific libraries:

* NumPy
* Pandas
* Matplotlib
* SciPy

For notebook execution, Jupyter and ipykernel are additionally required.

This repository is intended for scientific analysis, validation of correction schemes, and investigation of airmass-related retrieval effects in EM27/SUN observations. It is designed to be easily cloned and executed on different platforms, allowing consistent analysis workflows across multiple computers and operating systems.
