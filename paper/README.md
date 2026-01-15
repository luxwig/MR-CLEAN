# Paper Reproducibility (Scripts Used for the Manuscript)

This folder contains the MATLAB scripts and supporting functions used to generate the denoising results reported in the manuscript.

> Note: Some scripts may include study-specific configuration (e.g., file locations, scan-specific parameters). The long-term goal of this repository is to provide a more general-purpose toolbox interface; this `paper/` directory is kept to document and reproduce the results reported in the manuscript.

---

## Entry-point Scripts 
- **`asl.m`**: Denoising pipeline for the ASL Study (label/control data).
- **`cardiac.m`**: Denoising pipeline for the Cardiac Study.
- **`phantom.m`**: Denoising pipeline for the Phantom Study.

## Requirements

- MATLAB 
- BART Toolbox (v0.0.9, https://github.com/mrirecon/bart)
- GRAPPA Reconstruction Tools (https://github.com/mchiew/grappa-tools)
- mapVBVD (https://github.com/pehses/mapVBVD)