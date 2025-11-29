# A Geographically Reproducible Framework for Area of Habitat (AOH) Maps at National and Sub-National Scales

**Authors:** Wenxin Yang<sup>1</sup>, Lei Song<sup>2</sup>, Patrick R. Roehrdanz<sup>3</sup>, Nickolas McManus<sup>3</sup>, Dan Willett<sup>1</sup>, Elkin A. Noguera-Urbano<sup>4</sup>, Susana Rodríguez Buriticá<sup>4</sup>, Kevin Ramos<sup>1</sup>, Mary E. Blair<sup>5</sup>, Cristian A. Cruz-Rodríguez<sup>6</sup>, Carlos Jair Muñoz<sup>7</sup>, Amy E. Frazier<sup>1</sup>

<sup>1</sup>University of California, Santa Barbara, US,
<sup>2</sup>Rutgers University - New Brunswick, US,
<sup>3</sup>Conservation International, US,
<sup>4</sup>Instituto Humboldt, Colombia,
<sup>5</sup>American Museum of Natural History, US,
<sup>6</sup>University of Montreal, Canada,
<sup>7</sup>Lund University, Sweden.

<sup>*</sup>Corresponding author. Email: afrazier@ucsb.edu

**Recommended Citation:** To be added.

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
