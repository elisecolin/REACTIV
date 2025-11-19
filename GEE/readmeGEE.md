# 🌍 REACTIV – Google Earth Engine Scripts  
**Large-scale SAR change detection and visualization framework**  
*(Official GEE implementation for the paper: “REACTIV and Sentinel-1: Advances in Event-Based Alerts, Seasonal Adaptation, and Polarimetric Optimization for Large-Scale Change Detection”)*  

---

## 📖 Overview  

This folder contains the official **Google Earth Engine (GEE)** implementations of the four methodological extensions introduced in the paper.  
Each script is designed for direct execution within the **GEE Code Editor**, enabling large-scale visualization and rapid prototyping of Sentinel-1 time-series analyses.  

The REACTIV framework relies on the **Coefficient of Variation (CV)** as a simple yet robust descriptor of temporal dynamics in SAR backscatter.  
The extensions presented here enhance REACTIV’s interpretability and adaptability for different applications:  

| Extension | Purpose | Example Applications |
|------------|----------|----------------------|
| **Event-Based Detection (CV Ratio)** | Isolate newly occurring changes by comparing recent and historical temporal windows | Urban growth, ship detection, construction sites |
| **Frozen Background (Stable Reference)** | Build a temporally stable background by iteratively removing anomalies | Port activity, forest clearing, flood aftermath |
| **Seasonal Adaptation** | Restrict analysis to stable seasonal periods to filter cyclic natural variations | Agricultural monitoring, glacier dynamics |
| **Polarimetric Extension (MCV)** | Exploit multivariate CV across polarization channels to enhance weak-signal discrimination | RFI detection, glacial and cryospheric studies |

---

## ⚙️ Installation and Usage  

1. Open [**Google Earth Engine Code Editor**](https://code.earthengine.google.com/).  
2. Copy the desired script into your personal GEE workspace.  
3. Adjust the user parameters at the top of the script:  

   ```javascript
   var str2 = '2025-07-01';   // End date of the observation period
   var durationMonths = 6;    // Duration of the time series (months)
   var geometry = Map.getCenter(); // Area of interest

4. Run the script to generate REACTIV composites and detection layers.

5. Use the interactive chart panel (click on the map) to explore temporal amplitude profiles and cumulative CV behavior at any location.

Each script automatically selects the most frequently observed Sentinel-1 orbit and polarization configuration (VV/VH or single channel) to ensure geometric consistency.

🧩 Script List
Script	Description   
REACTIV_eventDetection	Detects new events using the Coefficient of Variation (CV) ratio between two temporal segments.  
REACTIV_frozenBackground	Builds a “frozen” background by progressively excluding outliers until statistical stability is reached.  
REACTIV_seasonal	Filters the time series to a given seasonal window (e.g., winter months) for cyclic change suppression.  
REACTIV_polarimetric	Implements the Multivariate Coefficient of Variation (MCV) for dual-polarization Sentinel-1 datasets.  
   
