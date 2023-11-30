# The-assessment-of-DG-blackout-risk-by-GCC

We provide the code and data relevant to the quality and nature of the data that is used to assess the system risk in the paper Assessment of Global Climate Change’s Impact on Electricity Distribution Service. 

Data
1) Climate prediction data from USGS NCCV database (US06037_MeanModel_english.csv)

   (https://www.usgs.gov/tools/national-climate-change-viewer-nccv) 

2) The fitted curves for hourly electricity load and temperature for 10 utilities.  (fitted_curve.mat)
   
   The local DG system data is sensitive information regarding the infrastructure security of the urban area and is subject to the non-disclosure agreement (NDA). We can not disclose all utility's hourly power consumption data. But we can disclose the fitted curves for electricity load and temperature for 10 utilities. They are then randomly distributed over the 33 nodes. These utilities were selected based on their availability and their representation of the overall load and temperature characteristics of the region. By providing the fitted curves for multiple utilities, we aimed to capture the variability in load and temperature patterns across the region. This information can be useful for operating power systems that are resilient to climate variability and change.

3) Distribution grid data (case33bw.m)
   
   In the paper, we used real distribution network data, but due to NDA restrictions we cannot make the data publicly available. So we demonstrated the calculation process using IEEE-33 as an example. We distributed fitted load temperature curves over the IEEE-33 case so that our code can be verified.
   
4) Hourly temperature data corresponding to the geo-location where the power consumption data is retrieved. (historical_temperature.mat)

5) Predicted load data between 1951-2100 for users at 33 nodes projected from the fitted load temperature curves and climate prediction data. (load_1951_2100_org.mat)

Code
1) The method to learn how consumer electricity behavior varies with temperature (prectict_load.m)
2) Risk assessment method (main.m)
3) The method of calculating the Jacobi matrix required to calculate the SOB (CalculateJacobian.m)
