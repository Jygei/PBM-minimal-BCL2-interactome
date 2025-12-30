[![DOI](https://zenodo.org/badge/973343500.svg)](https://doi.org/10.5281/zenodo.15288515)

# FLAME codes and analysis for *Stochasticity contributes to explaining minority and majority MOMP during apoptosis*

This repository contains the code and analysis files used to generate the results presented in the following publication:

> **Geiger, J., Klötzer, F., Pollak, N., Fullstone, G., Rehm, M. et al. (2025). _Stochasticity contributes to explaining minority and majority MOMP during apoptosis_. Cell Death & Disease, 16:893.**  
> DOI: https://doi.org/10.1038/s41419-025-08258-9

If you use this code in your work, please cite **both**:
1. the publication above, and  
2. this software (Zenodo DOI badge at the top).

---

## Citation

### Software
Geiger, J. (2025). *FLAME codes and analysis for Stochasticity contributes to explaining minority and majority MOMP during apoptosis*. Zenodo.  
DOI: https://doi.org/10.5281/zenodo.15288515

### Related publication
Geiger, J., Klötzer, F., Pollak, N., Fullstone, G., Rehm, M. et al. (2025).  
*Stochasticity contributes to explaining minority and majority MOMP during apoptosis*.  
Cell Death & Disease, 16:893.  
DOI: https://doi.org/10.1038/s41419-025-08258-9

---

# README for using the FLAME codes (see material and methods section):
We developed different particle-based models (whole-cell PBM and single-mito PBM) that can be used with FLAME v1.5, but the following files need to be changed:
- 0Creator.c
- AgentVariables.h
- functions.c
- getdata.c
- globals.h
- ReactionVariables.h
- ReadData.h
- WriteData.h
- XMLModelFile.xml

There was a problem with the random number generator (RNG) when the buffersize was changed. Increasing the buffersize in ```XMLModelFile.xml``` changed the results because the random seed is based on the declared buffersize variable.
Normally the RNG seed is defined by the maximum buffersize. The buffersize was statically set to 1045876.
This bug was fixed in FLAME v2, but since FaST was used to generate the code FLAME v1 was used. As a workaround the template of the header file and simulation.xslt (in ```..\FLAME-GPU\FLAMEGPU\templates```) have been changed.
Therefore, changes to the buffersize no longer have any effect on the RNG seed. The redefinition in the templates is marked/findable with "Jenny"

# README Python visualizer "pyvista_renderer" (see material and method section):
All visualizations of the particle-based model were created in Python 3.11. The integrated XML parser was used to parse the XML files of the particle-based model to create 3D rendered visualizations of the simulation results.

# README Python distribution analysis "heatmaps" (see materials and methods section):
All visualizations of the particle-based model were created in Python 3.11. The built-in XML parser was used to parse the XML files of the particle-based model to generate heatmaps showing the particle levels in terms of intensities.

---

Jenny Geiger  
Institute of Cell Biology and Immunology (IZI)  
University of Stuttgart, Germany



