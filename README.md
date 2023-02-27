These are datasets and codes for 'Statistical Modelling of Dengue Transmission Dynamics with Environmental Factors' and there are three files.

1: File 'Weather' contains dataset of weather factors for Southeast Asian and South American regions mentioned in the text and we have handled the missing values by replacing with the mean value. (air-pollution is not available for South American regions.)

2: File 'Dengue_number' contains dataset of dengue case count for all Southeast Asian and South American regions mentioned in the text and we have handled the missing values by replacing with the mean of the same year. (Originally, Singapore and South American regions have weekly data while Thailand has monthly data but we have converted the monthly data of Thailand to weekly data.)

3: File 'Code' contains R-code for the models and simulations in the project. 'all_step' is the code to illustrate the all-step estimation. 'P_value' is the code to calculate the p-values of each coefficient. You can first execute 'all_step' to import the packages.
