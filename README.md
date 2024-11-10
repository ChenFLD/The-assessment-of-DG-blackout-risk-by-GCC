# Assessment of Global Climate Change’s Impact on Electricity Distribution Service

This repository provides the MATLAB source code and a sample dataset for assessing power distribution system collapse risk under projected climate scenarios. The provided data and code simulate the assessment methods detailed in the accompanying paper.

## Project Structure

- `/src`: Contains all source code files necessary to run the assessment.
- `/data`: Contains a sample dataset for demonstration purposes, enabling users to test the code functionality without requiring external data.

## Source Code

The `/src` folder includes the following key scripts:

- **`predict_load.m`**: Predicts future electricity load based on temperature data. This script outputs a `.mat` file with predicted loads, which can be used in subsequent risk calculations.
- **`main.m`**: Calculates annual collapse risk for a power distribution network using the predicted load data and network configurations.
- **`CalculateJacobian.m`**: Computes the Jacobian matrix required to evaluate grid stability boundaries (SOB). This function supports input of different grid architectures.

## Sample Dataset

The `/data` folder includes files necessary to demonstrate the code’s capabilities:

- **`US06037_MeanModel_english.csv`**: Climate prediction data from the USGS NCCV database.
- **`fitted_curve.mat`**: Fitted curves representing hourly electricity load and temperature data for 10 utilities, used to simulate load variability.
- **`case33bw.m`**: IEEE-33 bus system configuration used in the demo. This can be replaced with other configurations for different network architectures.
- **`historical_temperature.mat`**: Historical temperature data for the region of interest.
- **`load_1951_2100_org.mat`**: Projected load data for each of the 33 nodes from 1951 to 2100, based on climate projections.

## System Requirements

- **MATLAB Version**: R2021a or later
- **Dependencies**:
  - [CVX](http://cvxr.com/cvx/) for optimization
  - [MATPOWER](https://matpower.org/) for power flow analysis

## Installation Guide

1. Install MATLAB and required dependencies:
   - Follow the [CVX installation guide](http://cvxr.com/cvx/doc/install.html).
   - Follow the [MATPOWER installation instructions](https://matpower.org/docs/MATPOWER-manual.pdf).
2. Clone this repository and add it to your MATLAB path:
   ```bash
   git clone [https://github.com/ChenFLD/The-assessment-of-DG-blackout-risk-by-GCC/tree/version2]
Estimated Installation Time: ~5 minutes
