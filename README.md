[![DOI](https://zenodo.org/badge/973343500.svg)](https://doi.org/10.5281/zenodo.15288515)

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
Jenny Geiger<br>
Institute of Cell Biology and Immunology (IZI) at the University of Stuttgart, Germany
