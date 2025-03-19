# REACTIV
## Rapid and EAsy Change detection in radar TIme-series by Variation coefficient

![Logo](https://raw.githubusercontent.com/elisekoeniguer/REACTIV/master/REACTIV.png)

### Project Overview
REACTIV was developed at **Onera** as part of the MEDUSA project in 2018. It is designed to visualize stacks of SAR images with change detection highlighting. 
At the **End of 2024**: Major updates developed as part of the **MUSIC** chaire include:
  - A **polarimetric version** for multivariate variation coefficient analysis, implemented in both **Python** and **JavaScript**.
  - A **seasonal version**, as well as **FrozenBackGround** and **NewEvent**, available exclusively on **Google Earth Engine (GEE)**.
    
REACTIV exists in multiple programming languages:
- **MATLAB** ((no longer maintained)
- **Python 3** (for standalone and scientific analysis)
- **JavaScript GEE** (for the Google Earth Engine GEE platform)
- **EOBROWSER**: [Sentinel-1 REACTIV script](https://custom-scripts.sentinel-hub.com/sentinel-1/reactiv/)


## Installation
For GEE: clone the Git repository by running the following command in a terminal:
```bash
git clone https://earthengine.googlesource.com/users/elisecolinaa/REACTIV
```
For Python: a Python notebook is available

## VERSIONS
The repository contains multiple versions of the REACTIV algorithm:

### 1. CLASSIC
**Description:**
The classic version of REACTIV, based on the coefficient of variation for change detection in SAR time-series of urban areas.

**Reference:**
Colin, E., & Nicolas, J. M. (2020). Change detection based on the coefficient of variation in SAR time-series of urban areas. *Remote Sensing, 12*(13), 2089. [DOI:10.3390/rs12132089](https://doi.org/10.3390/rs12132089)

---
### 2. FrozenBackGround *(GEE only)*
**Description:**
Reimplementation of the method for detecting ephemeral objects in SAR time-series using frozen background-based change detection.

**Reference:**
Taillade, T., Thirion-Lefevre, L., & Guinvarc’h, R. (2020). Detecting ephemeral objects in SAR time-series using frozen background-based change detection. *Remote Sensing, 12*(11), 1720. [DOI:10.3390/rs12111720](https://doi.org/10.3390/rs12111720)

---
### 3. NewEvent *(GEE only)*
**Description:**
This version detects changes occurring after a given time-series. Useful for monitoring the appearance of new events following a reference period.

---
### 4. POLARIMETRY *(Python & GEE)*
**Description:**
Implementation of change detection based on the bounds of multivariate coefficients of variation, calculated in their literal formulation.

**Reference:**
Colin, E., & Ossikovski, R. (2024). Towards a Unified Formalism of Multivariate coefficients of Variation: Application to analyzing polarimetric speckle time series. *Journal of the Indian Society of Remote Sensing, 52*(12), 2625-2636.

---
### 5. Seasons *(GEE only)*
**Description:**
A version that computes time-series by restricting them to a given season.
- Useful for the cryosphere to exclude snowfall periods and changes in the snow cover.
- Relevant for vegetation to eliminate seasonal changes that could distort the detection of anthropogenic or environmental changes.

## Usage
Each file can be executed independently based on user requirements. Specific dependencies and execution parameters are described in each source file.

## Citation
If you use this code in your research projects, please cite the following paper:
**Colin, E., & Nicolas, J. M. (2020). Change detection based on the coefficient of variation in SAR time-series of urban areas. Remote Sensing, 12(13), 2089.**

For more details, read this article: [How to visualize changes in a radar timeline](https://medium.com/@elisecolin/how-to-visualize-changes-in-a-radar-timeline-fb79ef526c1e)

Do not hesitate to contact the first author Elise Colin for further inquiries.

## Contribution
Contributions are welcome! Please submit a pull request or open an issue for any suggestions or improvements.

## License
This project is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)**.  
You are free to share and adapt this work for **non-commercial purposes**, provided that you give appropriate credit.

For more details, see the full license here: [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

