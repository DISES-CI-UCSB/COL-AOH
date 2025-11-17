# A Geographically Replicable Framework for Area of Habitat (AOH) Maps at National and Sub-National Scales

**Authors:** Wenxin Yang, Amy E. Frazier, Lei Song, Patrick R. Roehrdanz, Nickolas McManus, Dan Willett, Elkin A. Noguera-Urbano, Susana Rodríguez Buriticá

## Overview

This repository provides a workflow for generating Area of Habitat (AOH) maps for species using occurrence data, habitat preferences, and land cover classifications. The framework is designed to be geographically replicable and works with multiple taxonomic groups (amphibians, birds, mammals, reptiles).

## Workflow

Scripts should be run in numerical order:

1. **`01_prepData.R`** - Prepares species occurrence data and matches with land cover classifications
2. **`02_randomCI.R`** - Computes translation matrices using random sampling to convert IUCN habitat preferences to land cover classes
3. **`03_createHabitatLayer.R`** - Creates habitat layers from land cover data
4. **`04_createSpeciesAOH.R`** - Generates AOH maps (with uncertainty) for each species using range maps, habitat preferences, and elevation data
5. **`05_validation.R`** - Validates AOH maps against occurrence data
6. **`06_lv1Matrix.R`** - Creates level-1 classification matrices (optional)

## Additional Scripts

- **`funcs.R`** - Helper functions used across scripts
- **`info.R`** - Habitat and land cover classification metadata
- **`SensitivityViz.R`** - Visualization scripts for sensitivity analysis
- **`combine_csv.py`** - Python utility for combining CSV files
- **`migrate_biomodelos.py`** - Python script for migrating Biomodelos data

## Requirements

R packages: `sf`, `terra`, `here`, `rgrass`, `dplyr`, `data.table`, `rredlist`, `tidyr`, `Matrix`, `tidyverse`, `pROC`, `kableExtra`, `ggplot2`, `reshape2`, `parallel`

Python packages: Standard libraries (for utility scripts)

## Usage

1. Ensure all required data files are in the `data/` directory
2. Update file paths in scripts as needed
3. Run scripts sequentially from `01_prepData.R` through `06_lv1Matrix.R`
