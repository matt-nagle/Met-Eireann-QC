import geopandas
import pandas as pd
import netCDF4
import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable

from Live_Station_Data_Importation import rain_wow, rain_wow_hourly_obs, temp_wow, temp_official, rain_official


###########################################################################################################
###########################################################################################################
#####################################    STATIC DATA IMPORTATION    #######################################
###########################################################################################################
###########################################################################################################


### MUST BE EDITED! ###
# Change file path to wherever the static data (boundaries, elevation map) are stored

base_path = r"/home/matt/Documents/Met_Eireann/Data/"

### MUST BE EDITED! ###



##### Import Boundaries Data #####

# Import and define the boundaries for the counties of the island of Ireland
path_to_ROI_boundaries_data = base_path + r"Boundaries/Republic of Ireland Counties - shp/7829fd91-e0c3-4246-8fb5-100f3b62f8272020329-1-8mtzfa.iikfy.shp"
path_to_NI_boundaries_data = base_path + r"Boundaries/Northern Ireland Counties - shp/OSNI_Open_Data_-_50K_Boundaries_-_NI_Counties.shp"

ROI_counties = geopandas.read_file(path_to_ROI_boundaries_data)
NI_counties = geopandas.read_file(path_to_NI_boundaries_data)

# Note that we will use the Irish Grid projection so that we can get accurate distances
ROI_counties = ROI_counties.to_crs({'init': 'epsg:29902'}) # Convert maps to the same coordinate reference system
NI_counties = NI_counties.to_crs({'init': 'epsg:29902'})  # Convert maps to the same coordinate reference system





##### Import Elevation Data #####

path_to_elevation_data = base_path + r"Maps/section_with_edge_to_grib_final02.grib1.nc"

# Import the elevation data with a 30 metre resolution
raw_elevations_30m = netCDF4.Dataset(path_to_elevation_data)

# Extract the latitude, longitude and altitude from the raw data
elevation_data_lat = raw_elevations_30m['latitude'][:]
elevation_data_long = raw_elevations_30m['longitude'][:]
elevation_data_alt = raw_elevations_30m['altitude'][:]

# Define latitude and longitude as numpy arrays
elevation_data_lat_array = np.asarray(elevation_data_lat)
elevation_data_long_array = np.asarray(elevation_data_long)







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






##### Define a function to add the elevation data to a geo dataframe #####

def add_elevation(gdf_of_interest, ref_lat = elevation_data_lat_array, ref_long = elevation_data_long_array):
    """
    This function adds an elevation data column to a geo dataframe
    
    Inputs: 
    1. gdf_of_interest is the geodataframe that the Altitude data will be added to
    2. ref_lat is the array of reference latitudes, in our case the lat coordinates from the 30m elev data
    3. ref_lat is the array of reference longitudes, in our case the long coordinates from the 30m elev data
    
    Outputs:
    1. None, the Altitude column is added to the gdf_of_interest
    """  
    station_altitudes = []
    for index in gdf_of_interest.index:
        closest_lat_index = find_nearest_point_1D(ref_lat, gdf_of_interest.loc[index, "Latitude"], 
                                                  print_nearest_val=False)
        closest_long_index = find_nearest_point_1D(ref_long, gdf_of_interest.loc[index, "Longitude"], 
                                                   print_nearest_val=False)
        
        station_altitudes.append(elevation_data_alt[0, 0, closest_lat_index, closest_long_index])
        
        
    
    gdf_of_interest["Altitude"] = station_altitudes
    
    return







##### Define a function to isolate the data from the date, time and data of interest #####

def isolate_data_of_interest(day_of_interest, month_of_interest, year_of_interest, time_of_interest,
                             type_of_data = "Temperature", add_elevation_bool = True, 
                             remove_missing_val = True, cols_to_remove_missing_val = None):
    """
    This function isolates the data of interest.
    
    Inputs: 
    1. day_of_interest must be "DD"
    2. month_of_interest must be "MM"
    3. year_of_interest must be "YYYY"
    4. time_of_interest must be of the 24 hour format "HH" and will isolate data from the times between HH and HH+1
    5. type_of_data must be either "Rainfall" or "Temperature"
    6. add_elevation must be boolean True/False which dictates whether to add the elevation data or not
    7. remove_missing_val must be a boolean True/False which dictates whether to remove rows with missing values
    8. cols_to_remove_missing_val must be a list of columns to consider when looking for missing values,
       by default it will remove rows with ANY missing values however if you only care about Rainfall Accumulation
       then you can choose to only remove rows with missing values in that to preserve the maximum amount of data
    
    Outputs:
    1. geodataframe of wow data of interest
    2. geodataframe of official data of interest
    3. geodataframe of combined data of interest
    """  

    if(type_of_data == "Rainfall"):
        
        ##### WOW Data #####
        
        # Isolate the data at that date and time
        rain_wow_time_of_int = rain_wow.loc[(rain_wow["Day"] == day_of_interest) & 
                                                      (rain_wow["Month"] == month_of_interest) &
                                                      (rain_wow["Year"] == year_of_interest) & 
                                                      (rain_wow["Time"].str[:2] == time_of_interest),].copy()

        
        # In the rainfall accumulation case we want to keep only the last observation 
        # before the hour of interest is finished
        rain_wow_time_of_int_last_obs = rain_wow_time_of_int.drop_duplicates(subset="Site Id", 
                                                                                       keep="last").copy()
        
        
        # merge with dataset containing the hourly rainfall accumulation
        merged_dataset = pd.merge(rain_wow_time_of_int_last_obs, rain_wow_hourly_obs[["Id", "Rainfall Accumulation Hourly"]], on="Id")
        # Make sure there are no repeated indexes which would cause issues
        merged_dataset.reset_index(inplace=True) 

        # Convert data of interest into a GeoDataFrame for plotting
        gdf_wow = geopandas.GeoDataFrame(merged_dataset, 
                                         geometry=geopandas.points_from_xy(
                                             merged_dataset.Longitude, 
                                             merged_dataset.Latitude))
        
        gdf_wow.crs = {"init":"epsg:4326"} # initialise the dataframe to have a crs
        gdf_wow = gdf_wow.to_crs({'init': 'epsg:29902'}) # convert crs to Irish Grid Projection
        
        
        ##### Official Data #####
        
        # Isolate the data at that date and time
        rain_official_time_of_int = rain_official.loc[(rain_official["Day"] == day_of_interest) & 
                                                                (rain_official["Month"] == month_of_interest) &
                                                                (rain_official["Year"] == year_of_interest) & 
                                                                (rain_official["Time"].str[:2] == time_of_interest),].copy()
        
        # In the rainfall accumulation case we want to keep only the last observation 
        # before the hour of interest is finished
        rain_official_time_of_int_last_obs = rain_official_time_of_int.drop_duplicates(subset="stationid",
                                                                                                 keep="last").copy()

        # Convert data of interest into a GeoDataFrame for plotting
        gdf_official = geopandas.GeoDataFrame(rain_official_time_of_int_last_obs, 
                                              geometry=geopandas.points_from_xy(
                                                  rain_official_time_of_int_last_obs.Longitude, 
                                                  rain_official_time_of_int_last_obs.Latitude))
        
        gdf_official.crs = {"init":"epsg:4326"} # initialise the dataframe to have a crs
        gdf_official = gdf_official.to_crs({'init': 'epsg:29902'}) # convert crs to Irish Grid Projection
        
        
        ##### Combined Data #####
        
        # Combine data frames
        gdf_combined = gdf_wow.append(gdf_official)
        # Make sure there are no repeated indexes which would cause issues
        gdf_combined.reset_index(inplace=True) 
        
        
    elif type_of_data == "Temperature":
        
        ##### WOW Data #####
        
        # Isolate the data at that date and time
        temp_wow_time_of_int = temp_wow.loc[(temp_wow["Day"] == day_of_interest) & 
                                                      (temp_wow["Month"] == month_of_interest) &
                                                      (temp_wow["Year"] == year_of_interest) & 
                                                      (temp_wow["Time"].str[:2] == time_of_interest),].copy()
        
        # In the temperature case we want to keep only the last observation 
        # before the hour of interest is finished
        temp_wow_time_of_int_last_obs = temp_wow_time_of_int.drop_duplicates(subset="Site Id", 
                                                                                       keep="last").copy()

        # Convert data of interest into a GeoDataFrame for plotting
        gdf_wow = geopandas.GeoDataFrame(temp_wow_time_of_int_last_obs,
                                         geometry=geopandas.points_from_xy(
                                             temp_wow_time_of_int_last_obs.Longitude, 
                                             temp_wow_time_of_int_last_obs.Latitude))
        
        gdf_wow.crs = {"init":"epsg:4326"} # initialise the dataframe to have a crs
        gdf_wow = gdf_wow.to_crs({'init': 'epsg:29902'}) # convert crs to Irish Grid Projection
        
        
        ##### Official Data #####
        
        # Isolate the data at that date and time
        temp_official_time_of_int = temp_official.loc[(temp_official["Day"] == day_of_interest) & 
                                                                (temp_official["Month"] == month_of_interest) &
                                                                (temp_official["Year"] == year_of_interest) & 
                                                                (temp_official["Time"].str[:2] == time_of_interest),].copy()
        
        # In the temperature case we want to keep only the last observation 
        # before the hour of interest is finished
        temp_offical_time_of_int_last_obs = temp_official_time_of_int.drop_duplicates(subset="stationid",
                                                                                                keep="last").copy()

        # Convert data of interest into a GeoDataFrame for plotting
        gdf_official = geopandas.GeoDataFrame(temp_offical_time_of_int_last_obs,
                                              geometry=geopandas.points_from_xy(
                                                  temp_offical_time_of_int_last_obs.Longitude, 
                                                  temp_offical_time_of_int_last_obs.Latitude))
        
        gdf_official.crs = {"init":"epsg:4326"} # initialise the dataframe to have a crs
        gdf_official = gdf_official.to_crs({'init': 'epsg:29902'}) # convert crs to Irish Grid Projection
        
        
        
        ##### Combined Data #####
        
        # Combine data frames
        gdf_combined = gdf_wow.append(gdf_official)
        # Make sure there are no repeated indexes which would cause issues
        gdf_combined.reset_index(inplace=True) 
        
        
        
    ##### Add Elevation Data #####
    
    if(add_elevation_bool):
        add_elevation(gdf_of_interest = gdf_wow,
                      ref_lat = elevation_data_lat_array,
                      ref_long = elevation_data_long_array)
        add_elevation(gdf_of_interest = gdf_official,
                      ref_lat = elevation_data_lat_array,
                      ref_long = elevation_data_long_array)
        add_elevation(gdf_of_interest = gdf_combined,
                      ref_lat = elevation_data_lat_array,
                      ref_long = elevation_data_long_array)
        
    
    ##### Remove Missing Values #####
    
    if(remove_missing_val):
        gdf_wow.dropna(subset = cols_to_remove_missing_val, inplace = True)
        gdf_official.dropna(subset = cols_to_remove_missing_val, inplace = True)
        gdf_combined.dropna(subset = cols_to_remove_missing_val, inplace = True)
        
        
        
    return gdf_wow, gdf_official, gdf_combined 









    
##### Define a function that plots the data of interest #####


def plot_wow_data(gdf_of_interest, type_of_plot = "Air Temperature", 
                  buffer_val = 0, flags = None, end_of_title = None):
    """
    This function plots the data of interest
    
    Inputs: 
    1. gdf_of_interest containing the data of interest
    2. type_of_plot must be either "Rainfall Accumulation" or "Air Temperature"
    3. buffer_val (in metres) will plot a disk of radius buffer_val around each station 
    4. flags is the output from a QC check, any flags of 1 will plot a red dot on the station
    5. end_of_title is an optional argument that allows you to add additional text to the title
    
    Outputs:
    1. Desired plot
    """
    
    
    ##### Plot county borders for the Island of Ireland #####
    ax = ROI_counties["geometry"].plot(figsize=(10,15), edgecolors="grey", color="w")
    divider = make_axes_locatable(ax) # for vertically aligning the plot and the legend
    NI_counties["geometry"].plot(ax=ax, edgecolors="grey", color="w")

    if(end_of_title == None):
        ax.set_title(gdf_of_interest["Day"].iloc[0] + "-" + gdf_of_interest["Month"].iloc[0] + "-" + 
                     gdf_of_interest["Year"].iloc[0] + " " + gdf_of_interest["Time"].iloc[0][:2] + "am - " + 
                     type_of_plot)
    else:
        ax.set_title(gdf_of_interest["Day"].iloc[0] + "-" + gdf_of_interest["Month"].iloc[0] + "-" + 
                     gdf_of_interest["Year"].iloc[0] + " " + gdf_of_interest["Time"].iloc[0][:2] + "am - " + 
                     type_of_plot + " - " + end_of_title)        
    
    cax = divider.append_axes("right", size="5%", pad=0.1)
    
    
    
    ##### Plot Markers #####
    gdf_of_interest.plot(type_of_plot, ax=ax, legend = True, cax=cax, markersize=400, alpha = 0.5)

    # Annotate each marker with the Rainfall Accumulation to the nearest whole number
    for x, y, label in zip(gdf_of_interest.geometry.x, gdf_of_interest.geometry.y, gdf_of_interest[type_of_plot]):
        ax.annotate(round(label), xy=(x, y), verticalalignment='center', horizontalalignment='center', 
                    weight = "bold", fontsize = "large")
    

    ##### Plot of Buffer #####
    if(isinstance(buffer_val, int)):
        if(buffer_val > 0):
            new_gdf_of_interest = gdf_of_interest.copy()
            new_gdf_of_interest['geometry'] = new_gdf_of_interest['geometry'].buffer(buffer_val)

            new_gdf_of_interest["geometry"].plot(ax=ax, alpha = 0.2)
    elif(isinstance(buffer_val, list)):
        new_gdf_of_interest = gdf_of_interest.copy()
        new_gdf_of_interest['buffer'] = buffer_val
        new_gdf_of_interest['geometry'] = new_gdf_of_interest['geometry'].buffer(distance = new_gdf_of_interest['buffer'])
        
        new_gdf_of_interest = new_gdf_of_interest.loc[new_gdf_of_interest["buffer"] != 0,:]

        new_gdf_of_interest['geometry'].plot(ax=ax, alpha = 0.2)
    else:
        print("Error: Type of buffer_val must be integer or list!")
    
            
    ##### Flag Stations #####
    if flags is not None:
        flag_gdf_of_interest = gdf_of_interest.copy()
        flag_gdf_of_interest['flag'] = flags
        
        # Flag bad values red
        if(1 in flags):
            print("red: Stations that have been flagged as bad")
            bad_values = flag_gdf_of_interest.loc[flag_gdf_of_interest["flag"] == 1,].copy()
            bad_values["geometry"].plot(ax=ax, markersize=400, color="red")

        # Flag isolated inner values purple
        if(11 in flags):
            print("orange: Stations that are inner circle isolated (< 2 neighbours)")
            isolated_inner = flag_gdf_of_interest.loc[flag_gdf_of_interest["flag"] == 11,].copy()
            isolated_inner["geometry"].plot(ax=ax, markersize=400, color="orange")
        
        # Flag isolated outer values magenta
        if(12 in flags):
            print("pink: Stations that are outer circle isolated (< num_min_outer neighbours)")
            isolated_outer = flag_gdf_of_interest.loc[flag_gdf_of_interest["flag"] == 12,].copy()
            isolated_outer["geometry"].plot(ax=ax, markersize=400, color="hotpink")
            
        if(-999 in flags):
            print("grey: Stations that have not been checked")
            isolated_outer = flag_gdf_of_interest.loc[flag_gdf_of_interest["flag"] == -999,].copy()
            isolated_outer["geometry"].plot(ax=ax, markersize=400, color="dimgrey")
