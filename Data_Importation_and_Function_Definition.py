import geopandas
import pandas as pd
import netCDF4
import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable




###########################################################################################################
###########################################################################################################
#########################################    DATA IMPORTATION    ##########################################
###########################################################################################################
###########################################################################################################


##### Import Boundaries Data #####

base_path = r"/home/matt/Documents/Met_Eireann/Data/"

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





##### Import WOW Data #####

# Import Raw Wow data for the month of May
path_to_raw_wow_data = base_path + r"WOW_May_2021.csv"

raw_wow = pd.read_csv(path_to_raw_wow_data, dtype = {"Ground State":"string"})

# Create Day, Month, Year and Time columns
tmp = raw_wow["Report Date / Time"].str.split("-", expand=True)
tmp2 = tmp.iloc[:,2].str.split(" ", expand=True)

raw_wow["Day"] = tmp2.iloc[:,0]
raw_wow["Month"] = tmp.iloc[:,1]
raw_wow["Year"] = tmp.iloc[:,0]
raw_wow["Time"] = tmp2.iloc[:,1]

raw_wow.loc[raw_wow["Time"].isna(), "Time"] = "00:00:00" # Replace "missing" values with midnight

raw_wow["Hour"] = raw_wow["Time"].str[:2]

# Define a dataframe for the Rain Accumulation data
rain_wow = raw_wow[['Id', 'Site Id', 'Longitude', 'Latitude', 'Report Date / Time',
                    'Day', 'Month', 'Year', 'Time', 'Hour', 'Rainfall Accumulation', 'Rainfall Rate']].copy()


# Isolate the last observation every hour to reverse engineer the Hourly Rainfall Accumulation
rain_wow_hourly_obs = rain_wow.drop_duplicates(subset = ["Site Id", "Hour", "Day", "Month", "Year"], keep = "last").copy()

hourly_rainfall_accumulation = rain_wow_hourly_obs.groupby(["Site Id", "Day", "Month", "Year"])["Rainfall Accumulation"].diff().fillna(0)
rain_wow_hourly_obs["Rainfall Accumulation Hourly"] = hourly_rainfall_accumulation

# Define a dataframe for the Air Temperature data
temp_wow = raw_wow[['Id', 'Site Id', 'Longitude', 'Latitude', 'Report Date / Time',
                    'Day', 'Month', 'Year', 'Time', 'Hour', 'Air Temperature']].copy()






##### Import Official Data #####

##### Temperature Data #####

# Import Raw Official data for the month of May
path_to_raw_official_data = base_path + r"All stations May 2021 Full (updated lat long).csv"
raw_official = pd.read_csv(path_to_raw_official_data)


# Create Day, Month, Year and Time columns
tmp = raw_official["datein"].str.split("-", expand=True)
tmp2 = tmp.iloc[:,2].str.split("\xa0", expand=True)

raw_official["Day"] = tmp.iloc[:,0]
raw_official["Month"] = tmp.iloc[:,1]
raw_official["Year"] = tmp2.iloc[:,0]
raw_official["Time"] = tmp2.iloc[:,1]
raw_official["Hour"] = raw_official["Time"].str[:2]


# Define a dataframe for the Temperature data
temp_official = raw_official[['filetag', 'stationid', 'long', 'lat', 'datein',
                              'Day', 'Month', 'Year', 'Time', 'Hour', 't1dry', 't2dry']].copy()
# Replace values of -99.0 with missing values to be removed
temp_official.loc[temp_official["t1dry"] == -99.0, "t1dry"] = np.nan

temp_official.replace("May", "05", inplace=True, regex=True) # Replace month with number
temp_official.rename(columns={"t1dry":"Air Temperature", # Select t1dry as our temp obs
                              "long":"Longitude", 
                              "lat":"Latitude"}, inplace=True) 


##### Rainfall Data #####

# Import Raw Official rainfall data for the month of May
path_to_raw_official_rain_data = base_path + r"Rainfall stations May 2021.csv"
raw_official_rain = pd.read_csv(path_to_raw_official_rain_data)

# Create Day, Month, Year and Time columns
tmp = raw_official_rain["datein"].str.split("-", expand=True)
tmp2 = tmp.iloc[:,2].str.split("\xa0", expand=True)

raw_official_rain["Day"] = tmp.iloc[:,0]
raw_official_rain["Month"] = tmp.iloc[:,1]
raw_official_rain["Year"] = tmp2.iloc[:,0]
raw_official_rain["Time"] = tmp2.iloc[:,1]
raw_official_rain["Hour"] = raw_official_rain["Time"].str[:2]

# Define a dataframe for the Rain Accumulation data
rain_official = raw_official_rain[['filetag', 'stationid', 'long', 'lat', 'datein',
                                   'Day', 'Month', 'Year', 'Time', 'Hour',
                                   'totalpluvioaccrt_nrt','totalpluvioaccnrt', 'pluviohourrain']].copy()

# Define an additional column for the total rainfall accumulation rather than the accumulation over the last hour
rain_official["Rainfall Accumulation"] = rain_official.groupby(["stationid", "Day", "Month", "Year"])["totalpluvioaccnrt"].cumsum()

print("Note: As advised, the Official Rainfall Rate column is currently beign set as equal to the \nRainfall Accumulation Hourly column")
rain_official["Rainfall Rate"] = rain_official["totalpluvioaccnrt"]

rain_official.replace("May", "05", inplace=True, regex=True) # Replace month with number
rain_official.rename(columns={"totalpluvioaccnrt":"Rainfall Accumulation Hourly",
                              "long":"Longitude", 
                              "lat":"Latitude"}, inplace=True)








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
                  buffer_val = 0, flags = None):
    """
    This function plots the data of interest
    
    Inputs: 
    1. gdf_of_interest containing the data of interest
    2. type_of_plot must be either "Rainfall Accumulation" or "Air Temperature"
    3. buffer_val (in metres) will plot a disk of radius buffer_val around each station 
    4. flags is the output from a QC check, any flags of 1 will plot a red dot on the station
    
    Outputs:
    1. Desired plot
    """
    
    
    ##### Plot county borders for the Island of Ireland #####
    ax = ROI_counties["geometry"].plot(figsize=(10,15), edgecolors="grey", color="w")
    divider = make_axes_locatable(ax) # for vertically aligning the plot and the legend
    NI_counties["geometry"].plot(ax=ax, edgecolors="grey", color="w")

    ax.set_title(gdf_of_interest["Day"].iloc[0] + "-" + gdf_of_interest["Month"].iloc[0] + "-" + 
                 gdf_of_interest["Year"].iloc[0] + " " + gdf_of_interest["Time"].iloc[0][:2] + "am - " + 
                 type_of_plot)    
    
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
        flag_gdf_of_interest = flag_gdf_of_interest.loc[flag_gdf_of_interest["flag"] == 1,]

        flag_gdf_of_interest["geometry"].plot(ax=ax, markersize=400, color="red")

