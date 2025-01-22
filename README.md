# Codes used to generate figures of the article HERA: a high-resolution pan-European hydrological reanalysis (1951–2020)

This repository contains the R and python scripts used to assess performances of HERA and generate all figure in the article
Tilloy, A., Paprotny, D., Grimaldi, S., Gomes, G., Bianchi, A., Lange, S., Beck, H., and Feyen, L.: HERA: a high-resolution pan-European hydrological reanalysis (1950–2020), 
Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2024-41, in review, 2024.

## Table of Contents
- [Introduction](#introduction)
- [Scripts](#installation)
- [Dependencies](#dependencies)
- [License](#license)

## Introduction

The HERA high-resolution pan-European hydrological reanalysis (1951-2020) dataset [10.2905/a605a675-9444-4017-8b34-d66be5b18c95] is the result of a joint effort between the JRC and PIK to produce a long term hydrological reanalysis
with downscaled and bias-corrected climate reanalysis (ERA5-land) and dynamic socioeconomic inputs. It includes maps of climate variables (evaporation, evapotranspiration, precipitation, temperature),
dynamic socioeconomic inputs (land use, water demand, reservoir maps) required for hydrological modelling with LISFLOOD and river discharge with European extent at 1 arc minute (~1.5 km) grid resolution
and 6-hourly time step. TThis repository contains the scripts used to assess the peroformace of HERA against observed river discharge across Europe.

## Scripts 

The repository is divided into three main subfolders

1. [Validation]
2. [Usage]
3. [Supplement]

### Validation
- SpatialMatch.py: python script matching river gauges with HERA river pixels
- 
### Usage

[Provide guidelines on how to use the project. Include any necessary configuration, command-line options, or environment variables. Offer examples and use cases to help users understand the project's capabilities.]

## Dependencies

[List the external libraries, frameworks, or tools that your project depends on. Provide version numbers and links to their respective documentation or GitHub repositories.]

- [Dependency 1]
- [Dependency 2]
- [Dependency 3]

## Contributing

[Explain how others can contribute to the project. Include guidelines for reporting bugs, requesting features, and submitting pull requests. Mention any required code standards, coding practices, or testing procedures.]

## License

[State the license under which the project is released. Provide a link to the full license text, and explain any restrictions, obligations, or warranties associated with the license.]

---

[![GitHub Repository](https://img.shields.io/badge/GitHub-Repo%20Link-blue?style=flat&logo=github&logoColor=white)](https://github.com/Alowis/HERA)




# Title: Creation of the plots of the Usage note section of the article:  
# HERA: a high-resolution pan-European hydrological reanalysis (1950-2020)
# Author: Alois Tilloy - Joint Research Centre - Unit C6 
# Date: 2024 -02 -01 
# Description:
#   This script allows to generate plots that compare hydrological regimes 
#   of selected rivers at different points in time
# Note:
# Please note that the script requires transformed inputs from the HERA dataset
