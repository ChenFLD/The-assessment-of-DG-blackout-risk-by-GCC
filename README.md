# Assessment of Global Climate Change’s Impact on Electricity Distribution Service

This repository provides MATLAB source code and a sample dataset for assessing collapse risk in electricity distribution systems under projected climate scenarios. The demo uses the IEEE-33 bus system to illustrate the analysis. Users can also input their own network configurations and prediction data.

## Project Structure

- `/src`: Contains the source code files necessary for the assessment.
- `/data`: Contains the sample dataset used to demonstrate code functionality.
- `/demo`: Contains the necessary scripts and data files to reproduce the collapse risk assessment on the IEEE-33 bus system.

## Source Code

The `/src` folder includes the following main scripts:

- **`predict_load.m`**: Learns how temperature influences electricity load and generates future load scenarios based on input historical temperature, electricity load, and climate prediction data. Outputs a `.mat` file containing future load predictions.
- **`main.m`**: Computes annual collapse risk for a power distribution network. Calls the `CalculateJacobian.m` function to determine grid stability boundaries.
- **`CalculateJacobian.m`**: Calculates the Jacobian matrix, necessary for stability boundary assessments. This function can take different network architectures as inputs.

## Data

The `/data` folder includes files necessary to demonstrate code functionality:

- **`US06037_MeanModel_english.csv`**: Climate prediction data from the USGS NCCV database.
- **`fitted_curve.mat`**: Fitted load and temperature curves for 10 utilities, used to simulate variability across 33 nodes.
- **`case33bw.m`**: IEEE-33 bus system configuration used in the demo, replaceable with custom network configurations.
- **`historical_temperature.mat`**: Historical temperature data corresponding to the power consumption data’s location.
- **`load_1951_2100_org.mat`**: Projected load data for each node (1951–2100) based on climate projections.

## System Requirements

- **MATLAB Version**: R2021a or later
- **Dependencies**:
  - [CVX](http://cvxr.com/cvx/) for optimization
  - [MATPOWER](https://matpower.org/) for power flow analysis
- **Operating Systems**: Windows, macOS, and Linux

## Installation Guide

1. Install MATLAB and the required dependencies:
   - Follow the [CVX installation guide](http://cvxr.com/cvx/doc/install.html).
   - Follow the [MATPOWER installation guide](https://matpower.org/docs/MATPOWER-manual.pdf).
2. Clone this repository and add it to your MATLAB path:
   ```bash
   git clone https://github.com/YourUsername/ProjectName.git
   cd ProjectName
   ```
   Estimated Installation Time: ~10 minutes.

## Demo

The demo showcases the assessment method on the IEEE-33 bus configuration with provided data.

### Instructions to Run on Data

#### Run Risk Assessment:

Execute `main.m` to calculate the annual collapse risk based on the IEEE-33 bus configuration and projected load data. This script will call `CalculateJacobian.m` and use the provided network configuration (`case33bw.m`) and load data (`load_1951_2100_org.mat`).

```matlab
main
```

### Expected Output

- `main.m` outputs a set of collapse risk values for each year, based on the grid configuration.

### Expected Run Time

- **Risk Assessment**: ~10 minutes for IEEE-33 bus on a standard desktop.

## Instructions for Use

### Using Custom Network Architecture:

- Replace the default `case33bw.m` file in `/data` with your own network configuration.

### Using Custom Prediction Data:

- Replace `historical_temperature.mat` in the `/data` folder with your own historical electricity load data (due to NDA restrictions, we cannot disclose specific data; users can substitute their own data).
- Run `predict_load.m` in `/src` to perform step 1 and 2, generating `fitted_curve.mat` based on the historical temperature and load data.
- Replace `US06037_MeanModel_english.csv` with your own climate prediction data and use the generated `fitted_curve.mat` to produce future load predictions `load_1951_2100_org.mat`. This script will save the predictions in a `.mat` file.

```matlab
predict_load
```

- Rerun the `main.m` script following the same steps as in the demo to assess collapse risk with the new load data.
- Adjust any necessary parameters in `main.m` or `CalculateJacobian.m` to accommodate the new configuration.


## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.



