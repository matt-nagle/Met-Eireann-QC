import pandas as pd
import numpy as np


###########################################################################################################
###########################################################################################################
######################################    LIVE DATA IMPORTATION    ########################################
###########################################################################################################
###########################################################################################################


### MUST BE EDITED! ###
# Change file path to wherever the live data (wow station obs, official station obs) are stored
# Eventually this should be changed to directly interface with the live databases rather than
# pulling in csv

base_path = r"/home/matt/Documents/Met_Eireann/Data/"

### MUST BE EDITED! ###





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

