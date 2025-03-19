# REACTIV: Change Detection in SAR Time-Series

## Introduction
REACTIV is a collection of change detection algorithms applied to SAR (Synthetic Aperture Radar) time-series images. This project includes several implementations based on different scientific approaches to detect changes in various environments, such as urban areas, the cryosphere, and vegetation surfaces.

## Installation
Clone the Git repository by running the following command in a terminal:
```bash
git clone https://earthengine.googlesource.com/users/elisecolinaa/REACTIV
```

## Repository Contents
The repository contains five main files, each corresponding to a specific implementation of the REACTIV algorithm:

### 1. CLASSIC
**Description:**
The classic version of REACTIV, based on the coefficient of variation for change detection in SAR time-series of urban areas.

**Reference:**
Colin Koeniguer, E., & Nicolas, J. M. (2020). Change detection based on the coefficient of variation in SAR time-series of urban areas. *Remote Sensing, 12*(13), 2089. [DOI:10.3390/rs12132089](https://doi.org/10.3390/rs12132089)

---
### 2. FrozenBackGround
**Description:**
Reimplementation of the method for detecting ephemeral objects in SAR time-series using frozen background-based change detection.

**Reference:**
Taillade, T., Thirion-Lefevre, L., & Guinvarc’h, R. (2020). Detecting ephemeral objects in SAR time-series using frozen background-based change detection. *Remote Sensing, 12*(11), 1720. [DOI:10.3390/rs12111720](https://doi.org/10.3390/rs12111720)

---
### 3. NewEvent
**Description:**
This version detects changes occurring after a given time-series. Useful for monitoring the appearance of new events following a reference period.

---
### 4. POLARIMETRY
**Description:**
Implementation of change detection based on the bounds of multivariate coefficients of variation, calculated in their literal formulation.

**Reference:**
Colin, E., & Ossikovski, R. (2024). Towards a Unified Formalism of Multivariate coefficients of Variation: Application to analyzing polarimetric speckle time series. *Journal of the Indian Society of Remote Sensing, 52*(12), 2625-2636.

---
### 5. Seasons
**Description:**
A version that computes time-series by restricting them to a given season.
- Useful for the cryosphere to exclude snowfall periods and changes in the snow cover.
- Relevant for vegetation to eliminate seasonal changes that could distort the detection of anthropogenic or environmental changes.

## Usage
Each file can be executed independently based on user requirements. Specific dependencies and execution parameters are described in each source file.

## Contribution
Contributions are welcome! Please submit a pull request or open an issue for any suggestions or improvements.

# License
This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
You are free to share and adapt this work for *non-commercial* purposes, provided that you give appropriate credit.
For more details, see the full license here: CC BY-NC 4.0: https://creativecommons.org/licenses/by-nc/4.0/

