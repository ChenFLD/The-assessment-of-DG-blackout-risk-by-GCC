# The-assessment-of-DG-s-blackout-risk-by-GCC

We provide the code and data relevant to the quality and nature of the data that is used to assess the system risk in the paper Vulnerability and Adaptation of Power Grids for Global Warming Induced Electricity Demands. 

Data
1) Climate prediction data from USGS NCCV database (US06037_MeanModel_english.csv)

   (https://www.usgs.gov/tools/national-climate-change-viewer-nccv) 

2) The fitted curves for hourly electricity load and temperature for 10 utilities.  (fitted_curve.mat)
   
   The local DG system data is sensitive information regarding the infrastructure security of the urban area and is subject to the non-disclosure agreement (NDA). We can not disclose all utility's hourly power consumption data. But we can disclose the fitted curves for electricity load and temperature for 10 utilities. They are then randomly distributed over the 123 nodes. These utilities were selected based on their availability and their representation of the overall load and temperature characteristics of the region. By providing the fitted curves for multiple utilities, we aimed to capture the variability in load and temperature patterns across the region. This information can be useful for operating power systems that are resilient to climate variability and change.
   
3) Hourly temperature data corresponding to the geo-location where the power consumption data is retrieved. (historical_temperature.mat)

4) Distribution grid data (case123.mat)
   
   In the paper, we used real distribution network data, but due to NDA restrictions we cannot make the data publicly available. So we demonstrated the calculation process using IEEE-123 as an example. The data based on the 123 nodes was created by adding noise to real data in order to show similar patterns to the real data.

Code
1) The mothod to learn how consumer electricity behaviour varies with temperature (prectict_load.m)
2) Risk assessment method (main_1951_2099.m)
3) The method of calculating the Jacobi matrix required to calculate the SOB (CalculateJacobian.m)
