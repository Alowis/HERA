# Codes used to generate figures of the article HERA: a high-resolution pan-European hydrological reanalysis (1951–2020)

This repository contains the R and python scripts used to assess performances of HERA and generate all figure in the article
Tilloy, A., Paprotny, D., Grimaldi, S., Gomes, G., Bianchi, A., Lange, S., Beck, H., and Feyen, L.: HERA: a high-resolution pan-European hydrological reanalysis (1950–2020), 
Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2024-41, in review, 2024.

## Table of Contents

- [Introduction](#introduction)
- [Content](#Content)

## Introduction

The HERA high-resolution pan-European hydrological reanalysis (1951-2020) dataset [10.2905/a605a675-9444-4017-8b34-d66be5b18c95] is the result of a joint effort between the JRC and PIK to produce a long term hydrological reanalysis
with downscaled and bias-corrected climate reanalysis (ERA5-land) and dynamic socioeconomic inputs. It includes maps of climate variables (evaporation, evapotranspiration, precipitation, temperature),
dynamic socioeconomic inputs (land use, water demand, reservoir maps) required for hydrological modelling with LISFLOOD and river discharge with European extent at 1 arc minute (~1.5 km) grid resolution
and 6-hourly time step. TThis repository contains the scripts used to assess the peroformace of HERA against observed river discharge across Europe.

## Content 

The repository is divided into three main subfolders

1. [Validation]
2. [Usage]
3. [Supplement]

### Validation
#### python
- SpatialMatch.py: match river gauges with HERA river pixels.
- Correct_Qmatch.py: correction of gauge locations after automatic spatial matching.
- Valid_allY.py: generation of the performance metrics for each station (KGE,r, bias, ...).
- mHM_Match.py: match river gauges with mHM points.

#### R
- Domain_calib.R: Script to identify the extend of calibration in HERA (see article Supplement).
- Validation_Reanalysis.R: Main script used to generate plots in the article.

#### csv files
inputs for the spatial matchings

### Usage
#### R
- Data_usage.R: Script generating Figure 10 of the article.
  
### Supplenment
#### R
- Revision_HERA.R: Script comparing performances of HERA and mHM across Europe ( article Supplement) 

[![GitHub Repository](https://img.shields.io/badge/GitHub-Repo%20Link-blue?style=flat&logo=github&logoColor=white)](https://github.com/Alowis/HERA)

