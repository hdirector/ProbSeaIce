Analysis scripts for "Probabilistic Forecasting of the Arctic Sea Ice Edge with Contour Modeling" (Director, Raftery, Bitz 2020). To run many of these scripts, data that is archived elsewhere is required. What data is needed and how it should be formatted or processed is noted at the top of each file.

The files fitAndGen.R and fitAndGenLocal.R are the main scripts for contour modeling. The former is designed for use on a cluster and the latter is designed to be run locally. The script to compute the two-component mixture model is wght_EM.R. The file taskTable.R produces a .rda file giving all years, months, lead 
times, and training lengths to forecast. Other files produce reference forecasts (clim_ref.R, dPersis_ref.R). A number of files are also used to analyze results and produce figures and/or tables (brierScores.R, MCFTrainLengths.R, MCMCDiagnostics.R, mappingFigs.R, regionMaps.R, reliabilityDiagrams.R, visualizeForecasts.R).

The summaries folder contains intermediate data analysis results that can be used to produce tables and figures without downloading external data. See details at top of relevant scripts.
