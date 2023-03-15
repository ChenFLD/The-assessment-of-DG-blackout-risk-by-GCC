# The-assessment-of-DG-s-blackout-risk-by-GCC

We provide the code and data relevant to the quality and nature of the data that is used to assess the system risk in the paper Vulnerability and Adaptation of Power Grids for Global Warming Induced Electricity Demands. 

Data
1) Climate prediction data from USGS NCCV database (https://www.usgs.gov/tools/national-climate-change-viewer-nccv) (US06037_MeanModel_english.csv)
2) Hourly power consumption data of the consumers at every node in the DG.  (pnas_work_savepoint2.mat)
3) Hourly temperature data corresponding to the geo-location where the power consumption data is retrieved. (test_data.mat)

Code
1) The mothod to learn consumers’ power consumption behavior (poly_prectict.m)
2) Risk assessment method (main_45.m)
3) The method of calculating the Jacobi matrix required to calculate the SOB (CalculateJacobian.m)
