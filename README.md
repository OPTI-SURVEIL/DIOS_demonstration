# DIOS_demonstration
This repository contains code, data and intermediate results for the demonstration in DIOS framework paper.


## Code

* **DIOS_demonstration.R** generates the results and figures in the demonstration.

* **DIOS_demonstration_functions.R** contains functions used in DIOS_demonstration.R.

## Data

* **coord.sampled.csv** has x, y coordinates, risk factor X (Xvalue) and the log prevalence rates when ρ = 0.1 or 0.3 (logYValue0.1, logYValue0.3) for the 30 initial sites.

* **coord.unsampled.cvs** has x, y cooordinates and risk factor X (Xvalue) for the 70 unobserved sites.

## Intermediate results

These results can be used to generate the figures in the paper without running the optimization.

* **OFV 0.1.csv** has the ID of realizations (dataID), the ID of unobserved sites (siteID), and the OFV1 and OFV2 after including the corresponding unobserved site when ρ = 0.1.

* **OFV 0.3.csv** has the ID of realizations (dataID), the ID of unobserved sites (siteID), and the OFV1 and OFV2 after including the corresponding unobserved site when ρ = 0.3.

* **OFV 0.1 grid.csv** has the ID of realizations (dataID), the ID of each pixel (siteID), and the OFV1 and OFV2 after including the corresponding unobserved site when ρ = 0.1. This file need to be uncompressed with [7z](https://www.7-zip.org/) from **OFV 0.1 grid.7z**.

* **OFV 0.3 grid.csv** has the ID of realizations (dataID), the ID of each pixel (siteID), and the OFV1 and OFV2 after including the corresponding unobserved site when ρ = 0.3. This file need to be uncompressed with [7z](https://www.7-zip.org/) from **OFV 0.3 grid.7z**.

* **SA.OFV1.01.csv** has the IDs for the optimal 3-site combinations (BestID), the optimal OFV1 (BestVar), and the changes of OFV1 along the chain (criterion.result) when ρ = 0.1. Each row represents the result from one chain. 

* **SA.OFV1.03.csv** has the IDs for the optimal 3-site combinations (BestID), the optimal OFV1 (BestVar), and the changes of OFV1 along the chain (criterion.result) when ρ = 0.3. Each row represents the result from one chain. 

* **SA.OFV2.01.csv** has the IDs for the optimal 3-site combinations (BestID), the optimal OFV2 (BestVar), and the changes of OFV1 along the chain (criterion.result) when ρ = 0.1. Each row represents the result from one chain. 

* **SA.OFV2.03.csv** has the IDs for the optimal 3-site combinations (BestID), the optimal OFV2 (BestVar), and the changes of OFV1 along the chain (criterion.result) when ρ = 0.3. Each row represents the result from one chain. 

