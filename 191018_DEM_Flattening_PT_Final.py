"""
This script can be used to stretch a Digital Elevation Model
raster along the x-axis.

Author: Pratyush Tripathy
Email: pratkrt@gmail.com

Date: 18 October, 2019
"""



from pyrsgis import raster
import numpy as np
import math, sys, time
import shapefile as shp

#Set recurrsion depth limit to avoid error at a later stage
sys.setrecursionlimit(100000)

t1 = time.time()

# read the raster file
inFile = r"E:\DEM_Flattening\Sample_DEM.tif"
ds, elevationArr = raster.read(inFile, bands=1)
row, col = elevationArr.shape

# Defining tolerance is important to make sure that the switch between zero
# to actual elevation values at the edges don't mess up the things
tolerance = np.amin(elevationArr[elevationArr != np.amin(elevationArr)])

def cellCenterGenerator(inFile):
    global longitude_range, latitude_range, cell_size
    ds, elevationArr = raster.read(inFile)

    # Extract basic information from the input file
    row, col = elevationArr.shape
    ul_lon, xcell, xoffset, ul_lat, yoffset, ycell = ds.GetGeoTransform()
    ul_lon = round(ul_lon, 4)
    ul_lat = round(ul_lat, 4)
    cell_size = (xcell + abs(ycell)) / 2

    # Create Latitude and Longitude gradient array
    latitude_range = np.linspace(ul_lat, ul_lat + ycell*(row-1), row) # ycell is already negative
    longitude_range = np.linspace(ul_lon, ul_lon + xcell*(col-1), col)
    longitude_range += xcell / 2
    latitude_range += ycell / 2

    # Print the longitude and latitude range
    #print("Longitude ranges from %f to %f" % (longitude_range.min(), longitude_range.max()))
    #print("Latitude ranges from %f to %f" % (latitude_range.min(), latitude_range.max()))
    
    longitude, latitude = np.meshgrid(longitude_range, latitude_range)

    print("Longitude ranges from %f to %f" % (longitude.min(), longitude.max()))
    print("Latitude ranges from %f to %f" % (latitude.min(), latitude.max()))
   

    # loop through the dem and extract values
    for y in range(row):
        for x in range(col):
            yield(elevationArr[y, x], longitude[y, x], latitude[y, x])

# Define a function to calculate the actual distance betweek two points
def d_actual(point1, point2, cell_size, tolerance=tolerance):
    d_elev = abs(elevation_data[point1] - elevation_data[point2])
    #print(point1, elevation_data[point1], point2, elevation_data[point2])
    if d_elev < tolerance:
        return(abs(cell_size - math.sqrt(d_elev**2 + cell_size**2)))
    else:
        return(0)

# Define a recurssive function to update the new coordinates to the left and right
def update_loc_left(point, cell_size):
    lon, lat = point.split("_")
    lat = float(lat)
    left_lon = float(lon) - cell_size
    left_point = "%.4f_%.4f" % (left_lon, lat)
    if not float(left_lon) <= longitude_range.min():
        displacement = d_actual(point, left_point, cell_size)
        lon_location = np.where(longitude_range==float(left_lon))[0][0]
        for n in range(lon_location, 0, -1):
            temp_point = "%.4f_%.4f" % (longitude_range[n], lat)
            if n == n_left_lon - 1:
                updated_coordinates[temp_point] = "%.4f_%.4f" % (longitude_range[n]-displacement, lat)
                for m in range(lon_location-1, 0, -1):
                    temp_point_m = "%.4f_%.4f" % (longitude_range[m], lat)
                    updated_coordinates[temp_point_m] = temp_point_m
            else:
                temp_lon, temp_lat = updated_coordinates[temp_point].split("_")
                updated_coordinates[temp_point] = "%.4f_%.4f" % (float(temp_lon)-displacement, lat)
        update_loc_left(left_point, cell_size)
    else:
        pass

def update_loc_right(point, cell_size):
    lon, lat = point.split("_")
    lat = float(lat)
    right_lon = float(lon) + cell_size
    right_point = "%.4f_%.4f" % (right_lon, lat)
    if not float(right_lon) >= longitude_range.max():
        displacement = d_actual(point, right_point, cell_size)
        lon_location = np.where(longitude_range==float(right_lon))[0][0]
        for n in range(lon_location, col, 1):
            temp_point = "%.4f_%.4f" % (longitude_range[n], lat)
            # This can lead to a big mess here
            # If the number of columns in the raster is even, use n_right_lon + 1
            # If the number of columns in the raster is odd, use n_right_lon
            if n == n_right_lon + 1:
                updated_coordinates[temp_point] = "%.4f_%.4f" % (longitude_range[n]+displacement, lat)
                for m in range(lon_location+1, col, 1):
                    temp_point_m = "%.4f_%.4f" % (longitude_range[m], lat)
                    updated_coordinates[temp_point_m] = temp_point_m
            else:
                temp_lon, temp_lat = updated_coordinates[temp_point].split("_")
                updated_coordinates[temp_point] = "%.4f_%.4f" % (float(temp_lon)+displacement, lat)
        update_loc_right(right_point, cell_size)
    else:
        pass

# Get the center points of the raster and the elevation values
elevation_data = dict()
updated_coordinates = dict()
for elev, lon, lat in cellCenterGenerator(inFile):
    elevation_data["%.4f_%.4f" % (lon, lat)] = elev

# Get the center longitude value
col_center = math.floor(col/2)
center_longitude = longitude_range[col_center]

# Check if number of longitudes are even or odd
if col % 2 == 0:
    n_left_lon = col_center
    n_right_lon = col_center
    # To export reference raster
    outArray = elevationArr > 10000
    outArray[:, col_center] = 1
else:
    n_left_lon = col_center
    n_right_lon = col_center + 1
    # To export reference raster
    outArray = elevationArr > 10000
    outArray[:, col_center] = 1
    
print("\nBasic computations done. Flattening DEM now...")

# Loop through all the latitudes for center longitude
for center_latitude in latitude_range:
    center_point = "%.4f_%.4f" % (center_longitude, center_latitude)
    updated_coordinates[center_point] = center_point
    
    # Loop through all the points to the left
    update_loc_left(center_point, cell_size = cell_size)
    #print("Left side sorted.")
    
    # Loop through all the points to the right
    update_loc_right(center_point, cell_size = cell_size)
    #print("Right side sorted.")

print("\nUpdated all the left and right coordinates.")
print("\nTotal cells in input raster", len(elevation_data))    
print(len(updated_coordinates))

print("\nExporting shapefile...")

# Export the updated coordinates as point shapefile
outFile = inFile.replace(".tif", "_flattened_3.shp")
with shp.Writer(outFile) as w:
    fields = ['ID', 'oldLon', 'oldLat', 'Elevation']
    for f in fields:
        w.field(f, 'C')
    for n, record in enumerate(updated_coordinates):
        old_lon, old_lat = record.split("_")
        lon, lat = updated_coordinates[record].split("_")
        lon = float(lon)
        lat = float(lat)
        w.point(lon, lat)
        w.record(n, old_lon, old_lat, elevation_data[record])

# Export the center cell raster for reference
#raster.export(outArray, ds, filename=inFile.replace(".tif", "_centerReference.tif"), bands=1, dtype='int')

t2 = time.time()
print("Total time taken: %.2f minutes" % ((t2-t1)/60))

