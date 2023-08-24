# SimMCMC

This is a R code to perform simulation-based MCMC method to infer kinetic and delay parameters of a non-Markovian system.

## File Description
1. SimMCMC_functions_public.R
> A library of source functions that are needed to perform inference. 

2. SimMCMC_givendata_noise_model.R
> The main function performing inference using a given data. By default, it takes a data with the name "Example_data.csv."

3. SimMCMC_sim_noise_model_public.R
> The main function performing inference using simulated data. This functions does not require any given data as it generates simulation data by itself. 

4. Example_data.csv
> A tabular data contains time series representing the intensity of fluorescence. This matrix has 100 columns, and the length of each column 151. The 100 columns represent 100 different timeseries data. Each data is observed at time 0, 1, ..., 150. This example data can be used as an input of the function 'SimMCMC_givendata_noise_model.R'


