import titanlib
import numpy as np

# Import Formatted Data and functions from other Python script
from Data_Importation_and_Function_Definition import isolate_data_of_interest, plot_wow_data

# Import Custom Checks
from Custom_Checks import custom_range_check_climatology, plot_beta_vs_population_flagged




#===============================================================================================================#
#============================================ CHECKS TO PERFORM ================================================#
#===============================================================================================================#


##########################
##### Always Perform #####
##########################

perform_range_check_climatological = True
perform_metadata_check = True



##############################
##### Perform if Desired #####
##############################

# if there are specific upper/lower limits that we want to immediately 
# flag as suspicious we can use a range check however the climatological
# range check should already be getting rid of the outliers so this 
# might be superfluous
perform_range_check = True

# if there is an accepted 'neighbourhood' size that we think the spatial
# consistency checks would work on we can enable the isolation check and
# then perform whatever spatially dependent checks on the stations that
# pass
perform_isolation_check = True

if(perform_isolation_check):
    perform_buddy_event_check = False
    perform_buddy_check = True
    perform_fgt = True
    perform_sct = False
    perform_sct_resistant = True





#===============================================================================================================#
#============================================ OBTAIN LIVE DATA =================================================#
#===============================================================================================================#


###############################
##### Temperature Example #####
###############################


gdf_temp_wow, gdf_temp_official, gdf_temp_combined = isolate_data_of_interest(day_of_interest="20", 
                                                                              month_of_interest="05",
                                                                              year_of_interest="2021", 
                                                                              time_of_interest="10",
                                                                              type_of_data="Temperature", 
                                                                              add_elevation_bool=True,
                                                                              remove_missing_val=True, 
                                                                              cols_to_remove_missing_val=["Air Temperature"])


############################
##### Rainfall Example #####
############################


# By default let's discard any missing data from Rainfall Accumulation, Rainfall Accumulation Hourly
# or Rainfall Rate.

# Note if we were only interested in one of these values we could maximise our data points by only
# removing missing values from the column of interest
gdf_rain_wow, gdf_rain_official, gdf_rain_combined = isolate_data_of_interest(day_of_interest="05", 
                                                                              month_of_interest="05",
                                                                              year_of_interest="2021", 
                                                                              time_of_interest="10",
                                                                              type_of_data="Rainfall", 
                                                                              add_elevation_bool=True,
                                                                              remove_missing_val=True, 
                                                                              cols_to_remove_missing_val=["Rainfall Accumulation", 
                                                                                                          "Rainfall Accumulation Hourly", 
                                                                                                          "Rainfall Rate"])

###############################
##### Data Frame to Check #####
###############################

# Let's choose temperature as our example
gdf_to_check = gdf_temp_combined

data_to_check = "Air Temperature"







#===============================================================================================================#
#============================================= PERFORM CHECKS ==================================================#
#===============================================================================================================#



#================================= Perform Range Check Climatological ====================================#

if(perform_range_check_climatological):
    
    ##### Input Variables #####
    
    gdf_of_interest = gdf_to_check
    
    beta = 0.8
    type_of_check = "Temperature"
    add_columns = False
    
    ##### Input Variables #####

    flags = custom_range_check_climatology(gdf_of_interest, beta, type_of_check, add_columns)
    
    gdf_to_check["Climatological Flags"] = flags
    
    print("\n\n\n#=====================================================================#\n")
    print("#============= Flags from Climatological Range Check =================#\n")
    print("#=====================================================================#\n")
    print(flags)
    
    # Plot the percentage of the population flagged vs beta (to assist with tuning the correct value for beta)
    plot_beta_vs_population_flagged(gdf_to_check, min_beta = 0, max_beta = 4, step_beta = 0.1, type_of_data = type_of_check)





#========================================= Perform Metadata Check ========================================#

if(perform_metadata_check):
    
    ##### Input Variables #####
    
    long_np = np.array(gdf_to_check["Longitude"])
    lat_np = np.array(gdf_to_check["Latitude"])
    elev_np = np.array(gdf_to_check["Altitude"])

    points = titanlib.Points(lat_np, long_np, elev_np)

    check_lat = True
    check_lon = True
    check_elev = True
    check_laf = False
    
    ##### Input Variables #####

    flags = titanlib.metadata_check(points, check_lat, check_lon, check_elev, check_laf)
    
    gdf_to_check["Metadata Flags"] = flags
    
    print("\n\n\n#=====================================================================#\n")
    print("#=================== Flags from Metadata Check =======================#\n")
    print("#=====================================================================#\n")
    print(flags)





#========================================== Perform Range Check ==========================================#

if(perform_range_check):
    
    ##### Input Variables #####
    
    np_values = np.array(gdf_to_check[data_to_check]) # Values you want to check Rainfall/Temp

    # We have the choice between defining a global min and max
    # Or defining individual min and max values for each station
    # (potential use case: Checking measurement is within the ranges of the stations sensor)

    global_min = np.array([-5.0]) 
    global_max = np.array([35.0])
    
    ##### Input Variables #####

    flags = titanlib.range_check(np_values, global_min, global_max)
    
    gdf_to_check["Range Flags"] = flags
    
    print("\n\n\n#=====================================================================#\n")
    print("#===================== Flags from Range Check ========================#\n")
    print("#=====================================================================#\n")
    print(flags)





#======================================== Perform Isolation Check ========================================#

if(perform_isolation_check):
    
    ##### Input Variables #####
    
    long_np = np.array(gdf_to_check["Longitude"])
    lat_np = np.array(gdf_to_check["Latitude"])

    points = titanlib.Points(lat_np, long_np)
    
    radius = 8000 # Radius around each station to check
    num_min = 3 # Minimum number of neighbours
    
    ##### Input Variables #####

    flags = titanlib.isolation_check(points, num_min, radius)
    
    gdf_to_check["Isolation Flags"] = flags
    
    print("\n\n\n#=====================================================================#\n")
    print("#=================== Flags from Isolation Check ======================#\n")
    print("#=====================================================================#\n")
    print(flags)





#======================================= Perform Buddy Event Check =======================================#

if(perform_buddy_event_check):
    
    ##### Input Variables #####
    
    long_np = np.array(gdf_to_check["Longitude"])
    lat_np = np.array(gdf_to_check["Latitude"])

    points = titanlib.Points(lat_np, long_np) # Location of each station observation

    np_values = np.array(gdf_to_check[data_to_check]) # Values you want to check Rainfall/Temp

    radius = np.full(points.size(), 8000) # Radius of neighbourhood of points to consider
    num_min = np.full(points.size(), 3) # Minimum number of neighbours for a valid test


    event_threshold = 0.2 # Threshold above which the event is said to occur (i.e. rain happened)
    threshold = 0.25 # Fraction of other observations in the neighbourhood that must agree with the obs being tested
    max_elev_diff = 0 # Difference between elevations within the neighbourhood must not exceed this value
    elev_gradient = 0 # Value for linearly rescaling values for different elevations (e.g. for temp -0.0065 deg C per m)
    num_iterations = 5 # Number of iterations
    
    ##### Input Variables #####

    flags = titanlib.buddy_event_check(points, np_values, radius, num_min, event_threshold, 
                                       threshold, max_elev_diff, elev_gradient, num_iterations)
    
    gdf_to_check["Buddy Event Flags"] = flags
    
    print("\n\n\n#=====================================================================#\n")
    print("#================== Flags from Buddy Event Check =====================#\n")
    print("#=====================================================================#\n")
    print(flags)





#========================================== Perform Buddy Check ==========================================#

if(perform_buddy_check):
    
    ##### Input Variables #####
    
    long_np = np.array(gdf_to_check["Longitude"])
    lat_np = np.array(gdf_to_check["Latitude"])

    points = titanlib.Points(lat_np, long_np) # Location of each station observation

    np_values = np.array(gdf_to_check[data_to_check]) # Values you want to check Rainfall/Temp

    radius = np.full(points.size(), 10000) # Radius of neighbourhood of points to consider
    num_min = np.full(points.size(), 5) # Minimum number of neighbours for a valid test (else throw an error)


    # The buddy check flags an observations if |observation - avg of neighbours| normalized by the 
    # standard deviation in the circle is greater than a predefined threshold

    threshold = 3 # standard deviation threshold for flagging suspicious observations
    min_std = 1 # If the sd of observations in the neighbourhood is less than this value use this instead as the check
    max_elev_diff = 0 # Difference between elevations within the neighbourhood must not exceed this value
    elev_gradient = -0.0065 # Value for linearly rescaling values for different elevations (e.g. for temp -0.0065 deg C per m)
    num_iterations = 5 # Number of iterations
    
    ##### Input Variables #####

    flags = titanlib.buddy_check(points, np_values, radius, num_min, threshold, 
                                 max_elev_diff, elev_gradient, min_std, num_iterations)
    
    gdf_to_check["Buddy Flags"] = flags
    
    print("\n\n\n#=====================================================================#\n")
    print("#===================== Flags from Buddy Check ========================#\n")
    print("#=====================================================================#\n")
    print(flags)





#======================================= Perform First Guess Test ========================================#

if(perform_fgt):
    
    ##### Input Variables #####
    
    long_np = np.array(gdf_to_check["Longitude"])
    lat_np = np.array(gdf_to_check["Latitude"])
    elevs_np = np.array(gdf_to_check["Altitude"])

    points = titanlib.Points(lat_np, long_np, elevs_np)
    values = np.array(gdf_to_check[data_to_check], dtype=float)

    # if we want to check all observations we can use:
    # obs_to_check = np.ones(len(gdf_to_check))
    
    # else if we want to just check observations that passed the isolation check:
    failed_isolation_check = np.array(gdf_to_check["Isolation Flags"])
    obs_to_check = 1 - failed_isolation_check
    
    background_values = np.zeros(len(gdf_to_check))
    background_uncertainties = np.ones(len(gdf_to_check))
    background_elab_type = titanlib.MedianOuterCircle
    N = len(lat_np)
    num_min_outer = 3
    num_max_outer = 999 # Our goal is not to optimize speed here so let's include all observations
    inner_radius = 10000
    outer_radius = 20000
    num_iterations = 5
    num_min_prof = 0
    min_elev_diff = 100
    min_horizontal_scale = 250 
    max_horizontal_scale = 100000
    kth_closest_obs_horizontal_scale = 2
    tpos = np.ones(N) * 3 # Let's flag anything 3 standard deviations away from the background val
    tneg = np.ones(N) * 3 # Let's flag anything 3 standard deviations away from the background val


    values_mina = values.copy() - 20
    values_maxa = values.copy() + 20

    values_minv = values - 1
    values_maxv = values + 1
    debug = False
    
    ##### Input Variables #####

    flags = titanlib.fgt(points, values, obs_to_check, background_values, background_uncertainties, 
                         background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, 
                         num_iterations, num_min_prof, min_elev_diff, values_mina, values_maxa, values_minv, 
                         values_maxv, tpos, tneg, debug)
    
    gdf_to_check["FGT Flags"] = flags[0]
    
    print("\n\n\n#=====================================================================#\n")
    print("#=================== Flags from First Guess Test =====================#\n")
    print("#=====================================================================#\n")
    print(flags)





#=================================== Perform Spatial Consistency Test ====================================#

if(perform_sct):
    
    ##### Input Variables #####
    
    num_min = 3
    num_max = 100
    inner_radius = 10000
    outer_radius = 20000
    num_iterations = 5
    num_min_prof = 20
    min_elev_diff = 200
    min_horizontal_scale = 10000
    vertical_scale = 200

    long_np = np.array(gdf_to_check["Longitude"])
    lat_np = np.array(gdf_to_check["Latitude"])
    elevs_np = np.array(gdf_to_check["Altitude"])

    points = titanlib.Points(lat_np, long_np, elevs_np)
    values = np.array(gdf_to_check[data_to_check], dtype=float)

    pos = np.full(points.size(), 4)
    neg = np.full(points.size(), 4)
    eps2 = np.full(points.size(), 0.5)
    
    ##### Input Variables #####

    flags = titanlib.sct(points, values, num_min, num_max, inner_radius, outer_radius,
                        num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, vertical_scale,
                        pos, neg, eps2)
    
    gdf_to_check["SCT Flags"] = flags[0]
    
    print("\n\n\n#=====================================================================#\n")
    print("#=============== Flags from Spatial Consistency Test =================#\n")
    print("#=====================================================================#\n")
    print(flags)





#========================================= Perform SCT Resistant =========================================#

if(perform_sct_resistant):
    
    ##### Input Variables #####
    
    long_np = np.array(gdf_to_check["Longitude"])
    lat_np = np.array(gdf_to_check["Latitude"])
    elevs_np = np.array(gdf_to_check["Altitude"])

    points = titanlib.Points(lat_np, long_np, elevs_np)
    values = np.array(gdf_to_check[data_to_check], dtype=float)

    # if we want to check all observations we can use:
    # obs_to_check = np.ones(len(gdf_to_check))
    
    # else if we want to just check observations that passed the isolation check:
    failed_isolation_check = np.array(gdf_to_check["Isolation Flags"])
    obs_to_check = 1 - failed_isolation_check
    
    background_values = np.zeros(len(lat_np))
    background_elab_type = titanlib.MedianOuterCircle
    N = len(lat_np)
    num_min_outer = 3
    num_max_outer = 999 # Our goal is not to optimize speed here so let's include all observations
    inner_radius = 10000
    outer_radius = 20000
    num_iterations = 10
    num_min_prof = 1
    min_elev_diff = 100
    min_horizontal_scale = 250 
    max_horizontal_scale = 100000
    kth_closest_obs_horizontal_scale = 2
    vertical_scale = 200
    tpos = np.ones(N) * 3 # Let's flag anything 3 standard deviations away from the background val
    tneg = np.ones(N) * 3 # Let's flag anything 3 standard deviations away from the background val

    eps2 = np.ones(N) * 0.5

    values_mina = values.copy() - 20
    values_maxa = values.copy() + 20

    values_minv = values - 1
    values_maxv = values + 1
    debug = False
    
    ##### Input Variables #####

    flags = titanlib.sct_resistant(points, values, obs_to_check, background_values, background_elab_type, 
                                   num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, 
                                   num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, 
                                   kth_closest_obs_horizontal_scale, vertical_scale, values_mina, values_maxa, 
                                   values_minv, values_maxv, eps2, tpos, tneg, debug)
    
    gdf_to_check["SCT Res Flags"] = flags[0]
    
    print("\n\n\n#=====================================================================#\n")
    print("#===================== Flags from SCT Resistant ======================#\n")
    print("#=====================================================================#\n")
    print(flags)





#=================================== Check to see if any checks failed ===================================#

failed_check = np.zeros(len(gdf_to_check), dtype = int)

if("Climatological Flags" in gdf_to_check):
    failed_check = failed_check | np.array(gdf_to_check["Climatological Flags"])

if("Metadata Flags" in gdf_to_check):
    failed_check = failed_check | np.array(gdf_to_check["Metadata Flags"])
    
if("Range Flags" in gdf_to_check):
    failed_check = failed_check | np.array(gdf_to_check["Range Flags"])
    
# Note we are not including Isolation Check here as if a station passes all other checks but is
# isolated we will assume that the measurement is good.
    
if("Buddy Event Flags" in gdf_to_check):
    failed_check = failed_check | np.array(gdf_to_check["Buddy Event Flags"])
    
if("Buddy Flags" in gdf_to_check):
    failed_check = failed_check | np.array(gdf_to_check["Buddy Flags"])
    
if("FGT Flags" in gdf_to_check):
    failed_check = failed_check | (np.array(gdf_to_check["FGT Flags"]) == 1) | (np.array(gdf_to_check["FGT Flags"]) == 11) | (np.array(gdf_to_check["FGT Flags"]) == 12)

if("SCT Flags" in gdf_to_check):
    failed_check = failed_check | (np.array(gdf_to_check["SCT Flags"]) == 1) | (np.array(gdf_to_check["SCT Flags"]) == 11) | (np.array(gdf_to_check["SCT Flags"]) == 12)

if("SCT Res Flags" in gdf_to_check):
    failed_check = failed_check | (np.array(gdf_to_check["SCT Res Flags"]) == 1) | (np.array(gdf_to_check["SCT Res Flags"]) == 11) | (np.array(gdf_to_check["SCT Res Flags"]) == 12)


gdf_to_check["Failed Checks"] = failed_check





#===============================================================================================================#
#========================================== PLOT FLAGGED STATIONS ==============================================#
#===============================================================================================================#

##### Climatological Range Check Plot #####

if("Climatological Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["Climatological Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "Climatological Flags")


##### Metadata Check Plot #####

if("Metadata Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["Metadata Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "Metadata Flags")
    
    
##### Range Check Plot #####
    
if("Range Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["Range Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "Range Flags")
    

##### Isolation Check Plot #####

if("Isolation Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["Isolation Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "Isolation Flags")


##### Buddy Event Check Plot #####
    
if("Buddy Event Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["Buddy Event Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "Buddy Event Flags")


##### Buddy Check Plot #####
    
if("Buddy Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["Buddy Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "Buddy Flags")
    
    
##### First Guess Test Plot #####

if("FGT Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["FGT Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "FGT Flags")


##### Spatial Consistency Test Plot #####

if("SCT Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["SCT Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "SCT Flags")


##### SCT Resistant Check Plot #####

if("SCT Res Flags" in gdf_to_check):
    flags = np.array(gdf_to_check["SCT Res Flags"])
    plot_wow_data(gdf_to_check, type_of_plot = data_to_check, buffer_val = 0, 
                  flags = flags, end_of_title = "SCT Res Flags")

