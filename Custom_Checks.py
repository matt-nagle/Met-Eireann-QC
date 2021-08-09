import geopandas
import pandas as pd
import netCDF4
import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pyplot as plt




###########################################################################################################
###########################################################################################################
#########################################    DATA IMPORTATION    ##########################################
###########################################################################################################
###########################################################################################################


##### Import Normal Climate Ranges #####

base_path = r"/home/matt/Documents/Met_Eireann/Data/"


# Import and define the Normal Temperature and Rainfall Ranges
path_to_clim_norm_temp_mean = base_path + r"Climate Normals/Temperature/Tmean_grid_9120.csv"
path_to_clim_norm_temp_TN = base_path + r"Climate Normals/Temperature/TN_GR_normals9120.csv"
path_to_clim_norm_temp_TX = base_path + r"Climate Normals/Temperature/TX_GR_normals9120.csv"

path_to_clim_norm_rain_mean = base_path + r"Climate Normals/Rainfall/RainfallGrids_draft-Normals9120.csv"


# Import as CSV
clim_norm_temp_mean = pd.read_csv(path_to_clim_norm_temp_mean)
clim_norm_temp_TN = pd.read_csv(path_to_clim_norm_temp_TN)
clim_norm_temp_TX = pd.read_csv(path_to_clim_norm_temp_TX)

clim_norm_rain_mean = pd.read_csv(path_to_clim_norm_rain_mean)


# Convert into a GeoDataFrame
gdf_clim_norm_temp_mean = geopandas.GeoDataFrame(clim_norm_temp_mean.copy(), 
                                                 geometry=geopandas.points_from_xy(
                                                     clim_norm_temp_mean.east, 
                                                     clim_norm_temp_mean.north))
gdf_clim_norm_temp_TN = geopandas.GeoDataFrame(clim_norm_temp_TN.copy(),
                                               geometry=geopandas.points_from_xy(
                                                   clim_norm_temp_TN.east, 
                                                   clim_norm_temp_TN.north))
gdf_clim_norm_temp_TX = geopandas.GeoDataFrame(clim_norm_temp_TX.copy(), 
                                               geometry=geopandas.points_from_xy(
                                                   clim_norm_temp_TX.east, 
                                                   clim_norm_temp_TX.north))

gdf_clim_norm_rain_mean = geopandas.GeoDataFrame(clim_norm_rain_mean, 
                                                 geometry=geopandas.points_from_xy(
                                                     clim_norm_rain_mean.east, 
                                                     clim_norm_rain_mean.north))




##### Extract the x and y ticks for the temperature and rainfall grids #####

# Mean temperature values (note only for ROI at the moment)
clim_norm_mean_temp_y = gdf_clim_norm_temp_mean["north"] 
clim_norm_mean_temp_y_np = np.array(clim_norm_mean_temp_y)
clim_norm_mean_temp_y_np_unique = np.unique(clim_norm_mean_temp_y_np)

clim_norm_mean_temp_x = gdf_clim_norm_temp_mean["east"] 
clim_norm_mean_temp_x_np = np.array(clim_norm_mean_temp_x)
clim_norm_mean_temp_x_np_unique = np.unique(clim_norm_mean_temp_x_np)


# Mean minimum and maximum monthly temperature values
clim_norm_minmax_temp_y = gdf_clim_norm_temp_TN["north"] # Note grids for TN and TX are identical
clim_norm_minmax_temp_y_np = np.array(clim_norm_minmax_temp_y)
clim_norm_minmax_temp_y_np_unique = np.unique(clim_norm_minmax_temp_y_np)

clim_norm_minmax_temp_x = gdf_clim_norm_temp_TN["east"] # Note grids for TN and TX are identical
clim_norm_minmax_temp_x_np = np.array(clim_norm_minmax_temp_x)
clim_norm_minmax_temp_x_np_unique = np.unique(clim_norm_minmax_temp_x_np)


# Mean monthly rainfall values
clim_norm_rain_y = gdf_clim_norm_rain_mean["north"] 
clim_norm_rain_y_np = np.array(clim_norm_rain_y)
clim_norm_rain_y_np_unique = np.unique(clim_norm_rain_y_np)

clim_norm_rain_x = gdf_clim_norm_rain_mean["east"] 
clim_norm_rain_x_np = np.array(clim_norm_rain_x)
clim_norm_rain_x_np_unique = np.unique(clim_norm_rain_x_np)




###########################################################################################################
###########################################################################################################
#######################################    FUNCTION DEFINITION   ##########################################
###########################################################################################################
###########################################################################################################


##### Define a function to find the closest value in a 1D array #####

def find_nearest_point_1D(array_to_search, value, print_nearest_val = True):
    """
    This function finds the closest value in a 1D array
    
    Inputs: 
    1. array_to_search is the array to be searched
    2. value 
    3. print_nearest_val must be a boolean True/False 
    
    Outputs:
    1. index of the closest value in the array
    """  
    array_to_search = np.asarray(array_to_search)
    idx = (np.abs(array_to_search - value)).argmin()
    
    if(print_nearest_val):
        print("Nearest Value:",array_to_search[idx])
        
    return idx



##### Define a function to add the temperature min/max climate normal values #####

def add_temp_clim_normals(gdf_of_interest,
                          grid_of_minmax_temp_clim_norm_y = clim_norm_minmax_temp_y_np_unique, 
                          grid_of_minmax_temp_clim_norm_x = clim_norm_minmax_temp_x_np_unique, 
                          grid_of_mean_temp_clim_norm_y = clim_norm_mean_temp_y_np_unique, 
                          grid_of_mean_temp_clim_norm_x = clim_norm_mean_temp_x_np_unique):
    """
    This function adds three columns to the gdf_of_interest, one for the mean temperature, one for the mean
    monthly minimum temperature and the last one for the mean monthly maximum temperature at that station location 
    over the last 10 years.
    
    Inputs: 
    1. gdf_of_interest is the geodataframe that the temperature climate normals data will be added to
    2. grid_of_minmax_temp_clim_norm_y is the array of y coordinates, in our case the y ticks from the mean minimum/maximum temperature
    climate normals data
    3. grid_of_minmax_temp_clim_norm_x is the array of x coordinates, in our case the x ticks from the mean minimum/maximum temperature
    climate normals data
    4. grid_of_mean_temp_clim_norm_y is the array of y coordinates, in our case the y ticks from the mean monthly temperature
    climate normals data
    5. grid_of_mean_temp_clim_norm_x is the array of x coordinates, in our case the x ticks from the mean monthly temperature
    climate normals data
    
    Outputs:
    1. None, the columns are added to the gdf_of_interest
    """  
    mean_monthly_min = []
    mean_monthly_max = []
    mean_monthly_temp = []
    for index in gdf_of_interest.index:
        # Find the closest x and y grid points for the mean min/max temperatures
        closest_y_index = find_nearest_point_1D(grid_of_minmax_temp_clim_norm_y, 
                                                gdf_of_interest.geometry.y[index], 
                                                print_nearest_val=False)
        closest_x_index = find_nearest_point_1D(grid_of_minmax_temp_clim_norm_x, 
                                                gdf_of_interest.geometry.x[index], 
                                                print_nearest_val=False)
        
        # Find the closest x and y grid points for the monthly mean temperature
        closest_y_index_mean = find_nearest_point_1D(grid_of_mean_temp_clim_norm_y, 
                                                     gdf_of_interest.geometry.y[index], 
                                                     print_nearest_val=False)
        closest_x_index_mean = find_nearest_point_1D(grid_of_mean_temp_clim_norm_x, 
                                                     gdf_of_interest.geometry.x[index], 
                                                     print_nearest_val=False)
        
        
        
        # Find the month of interest and define the correct format for the different file formats
        month_of_interest = int(gdf_of_interest["Month"][0])
        min_month_of_int_format = "Tn_m" + str(month_of_interest)
        max_month_of_int_format = "Tx_m" + str(month_of_interest)
        mean_month_of_int_format = "Tm_m" + str(month_of_interest)
        
        # Append relevant climate normal data
        mean_monthly_min.append(gdf_clim_norm_temp_TN.loc[
            (gdf_clim_norm_temp_TN["east"] == grid_of_minmax_temp_clim_norm_x[closest_x_index]) &
            (gdf_clim_norm_temp_TN["north"] == grid_of_minmax_temp_clim_norm_y[closest_y_index]),
            min_month_of_int_format].values[0])
        
        mean_monthly_max.append(gdf_clim_norm_temp_TX.loc[
            (gdf_clim_norm_temp_TX["east"] == grid_of_minmax_temp_clim_norm_x[closest_x_index]) &
            (gdf_clim_norm_temp_TX["north"] == grid_of_minmax_temp_clim_norm_y[closest_y_index]),
            max_month_of_int_format].values[0])
        
        # NOTE: We currently do not have the NI data so assume that for NI stations
        # the mean value is in the middle of the min/max values.
        if (len(gdf_clim_norm_temp_mean.loc[
            (gdf_clim_norm_temp_mean["east"] == grid_of_mean_temp_clim_norm_x[closest_x_index_mean]) &
            (gdf_clim_norm_temp_mean["north"] == grid_of_mean_temp_clim_norm_y[closest_y_index_mean]),:]) == 0):
            mean_monthly_temp.append(np.mean([mean_monthly_min[-1], mean_monthly_max[-1]]))
        else:
            mean_monthly_temp.append(gdf_clim_norm_temp_mean.loc[
            (gdf_clim_norm_temp_mean["east"] == grid_of_mean_temp_clim_norm_x[closest_x_index_mean]) &
            (gdf_clim_norm_temp_mean["north"] == grid_of_mean_temp_clim_norm_y[closest_y_index_mean]),
            mean_month_of_int_format].values[0])
            
    
    gdf_of_interest["Mean Monthly Min Temp"] = mean_monthly_min
    gdf_of_interest["Mean Monthly Max Temp"] = mean_monthly_max
    gdf_of_interest["Mean Monthly Temp"] = mean_monthly_temp
    
    return




##### Define a function to add the rainfall climate normal values #####

def add_rain_clim_normals(gdf_of_interest,
                          grid_of_rain_clim_norm_y = clim_norm_rain_y_np_unique, 
                          grid_of_rain_clim_norm_x = clim_norm_rain_x_np_unique):
    """
    This function adds a column to the gdf_of_interest for the mean monthly rainfall at that station
    location over the last 10 years.
    
    Inputs: 
    1. gdf_of_interest is the geodataframe that the temperature climate normals data will be added to
    2. grid_of_rain_clim_norm_y is the array of y coordinates, in our case the y ticks from the rain climate normals data
    3. grid_of_rain_clim_norm_x is the array of x coordinates, in our case the x ticks from the rain climate normals data
    
    Outputs:
    1. None, the column is added to the gdf_of_interest
    """  
    mean_monthly_rainfall = []
    for index in gdf_of_interest.index:
        closest_y_index = find_nearest_point_1D(grid_of_rain_clim_norm_y, 
                                                gdf_of_interest.geometry.y[index], 
                                                print_nearest_val=False)
        closest_x_index = find_nearest_point_1D(grid_of_rain_clim_norm_x, 
                                                gdf_of_interest.geometry.x[index], 
                                                print_nearest_val=False)
        
        month_of_interest = int(gdf_of_interest["Month"][0])
        month_of_int_format = "m" + str(month_of_interest)
        
        mean_monthly_rainfall.append(gdf_clim_norm_rain_mean.loc[
            (gdf_clim_norm_rain_mean["east"] == grid_of_rain_clim_norm_x[closest_x_index]) &
            (gdf_clim_norm_rain_mean["north"] == grid_of_rain_clim_norm_y[closest_y_index]),
            month_of_int_format].values[0])
        
        
    
    gdf_of_interest["Mean Monthly Rainfall"] = mean_monthly_rainfall
    
    return







##### Define a function to perform the custom climatological range check #####

def custom_range_check_climatology(gdf_of_interest,
                                   beta = 1,
                                   type_of_check = "Temperature", 
                                   add_columns = True):
    """
    This function checks to see if the observed value aligns with the climate normals data at those station locations    
    
    Inputs: 
    1. gdf_of_interest is the geodataframe that contains the station observation data to be checked
    2. beta is the scaling factor. Larger values of beta make the check more lenient, smaller values of 
    beta make the check more strict
    3. type_of_check must be either "Temperature" or "Rainfall" which checks either the Air Temperature 
    observations or the Rainfall Accumulation
    4. add_columns must be a boolean True/False and it indicates whether or not the code should add the
    additional columns containing the adjusted Climate normals and the beta value to the gdf_of_interest
    
    Outputs:
    1. flags, a binary array indicating whether or not a station was flagged by the check
    
    """        
        
    if(type_of_check == "Temperature"):
        # check if the Climate Normals have been added to the gdf of interest
        # if they have not then add them
        if(("Mean Monthly Min Temp" not in gdf_of_interest) or 
        ("Mean Monthly Max Temp" not in gdf_of_interest) or 
        ("Mean Monthly Temp" not in gdf_of_interest)):
            add_temp_clim_normals(gdf_of_interest)
            

        # Pull out the observed values, mean values and the climate normals
        observed_values = np.array(gdf_of_interest["Air Temperature"])
        mean_values = np.array(gdf_of_interest["Mean Monthly Temp"])
        lower_limit = np.array(gdf_of_interest["Mean Monthly Min Temp"])
        upper_limit = np.array(gdf_of_interest["Mean Monthly Max Temp"])
        
        # Shift the climate normals by the mean so they are centred on zero
        lower_limit_shifted = lower_limit - mean_values
        upper_limit_shifted = upper_limit - mean_values
        
        # Rescale the allowable range by multiplying the shifted limits by beta
        lower_limit_rescaled = beta * lower_limit_shifted
        upper_limit_rescaled = beta * upper_limit_shifted
        
        # Shift the climate normals back to their original positions centered around the mean
        adjusted_lower_limit = lower_limit_rescaled + mean_values
        adjusted_upper_limit = upper_limit_rescaled + mean_values
        
        # Add the Adjusted limits to the gdf_of_interest
        if(add_columns):
            lower_lim_col_name = "Adjusted Lower Limit (beta =" + str(beta) + ")"
            upper_lim_col_name = "Adjusted Upper Limit (beta =" + str(beta) + ")"
            gdf_of_interest[lower_lim_col_name] = adjusted_lower_limit
            gdf_of_interest[upper_lim_col_name] = adjusted_upper_limit  
        
        
        # Check to see if the observed value is within the adjusted range
        flags = (observed_values < adjusted_lower_limit) | (observed_values > adjusted_upper_limit)
        # Convert from True/False to 1/0
        flags = flags.astype(int)
        
        
        
        
    # On advisement from Noel Fitzpatrick the rainfall check has been simplified to simply check if more than
    # a months average rainfall has fallen so far on the day of interest (note that the time of interest should
    # be close to midnight). I have also included the option to scale the months average rainfall by the factor
    # beta in case we want to make the check more/less strict.
    
    if(type_of_check == "Rainfall"):
        # check if the Climate Normals have been added to the gdf of interest
        # if they have not then add them
        if("Mean Monthly Rainfall" not in gdf_of_interest):
            add_rain_clim_normals(gdf_of_interest)
    
        
        # Pull out the observed values and the climate normals
        observed_values = np.array(gdf_of_interest["Rainfall Accumulation"])
        mean_monthly_rainfall = np.array(gdf_of_interest["Mean Monthly Rainfall"])
        
        # Scale the mean monthly rainfall by a factor of beta
        adjusted_mean_monthly_rainfall = beta * mean_monthly_rainfall
        
        # Add the Adjusted Monthly Rainfall to the gdf_of_interest
        if(add_columns):
            adjusted_mean_col_name = "Adjusted Monthly Rainfall (beta =" + str(beta) + ")"
            gdf_of_interest[adjusted_mean_col_name] = adjusted_mean_monthly_rainfall
        
        
        # Check to see if the observed value is greater than the scaled monthly rainfall
        flags = (observed_values > adjusted_mean_monthly_rainfall)
        # Convert from True/False to 1/0
        flags = flags.astype(int)
        
    
    return flags







##### Define a function to plot the percentage of stations flagged vs the value of beta #####

def plot_beta_vs_population_flagged(gdf_of_interest, min_beta = 0, max_beta = 4, step_beta = 0.1, type_of_data = "Temperature"):
    """
    This function plots the percentage of the population of stations flagged vs the value of beta.
    
    This should help with tuning the value for beta as according to the Eischeid paper (see below) 
    the correct value of beta should be the value when the slope of the plot is approximately zero

    Ref: https://doi.org/10.1175/1520-0450(1995)034%3C2787:TQCOLT%3E2.0.CO;2
    
    Inputs: 
    1. gdf_of_interest is the geodataframe that contains the station observation data to be checked
    2. min_beta is the smallest value of beta to use
    3. max_beta is the largest value of beta to use
    4. step_beta is the step size to use for plotting the values of beta
    5. type_of_data must be either "Temperature" or "Rainfall" which flags the type of data being used
    
    Outputs:
    1. None, a plot is created.
    
    """
    beta = min_beta

    num_steps = (max_beta - min_beta)/step_beta

    beta_list = np.linspace(min_beta, max_beta, int(np.ceil(num_steps) + 1))
    perc_pop_flagged = []

    for beta in beta_list:
        flags = custom_range_check_climatology(gdf_of_interest, beta = beta, type_of_check = type_of_data, add_columns = False)

        perc_pop_flagged.append(sum(flags)/len(flags))


    # Plot percentage of population flagged vs beta
    plt.plot(beta_list, perc_pop_flagged)
    plt.xlabel("Beta")
    plt.ylabel("Percentage of Population Flagged")
    
    return