# Wave Attenuation Toolbox (WATTE) version 2
# Developed by Center for Computation & Technology and Center for Coastal Resiliency (Currently, Coastal Ecosystem Design Studio) at Louisiana State University (LSU).
# WATTE version 1 was originally developed by M. Foster-Martinez, University of New Orleans: Karim Alizad, University of South Carolina: Scott C. Hagen, LSU (https://digitalcommons.lsu.edu/civil_engineering_data/1/)
# Developer: Jin Ikeda, Shu Gao, Christopher E. Kees, and Peter Bacopoulos
# Deepest thanks: This software is dedicated to in memory of Dr. Scott C. Hagen. We are truly blessed to have been one of your pupils. We will spread your idea further.
#
# WATTE version 2 is an open source-based toolbox using Python 3 and QGIS (Verified at version Python 3.10 and QGIS 3.22). It estimates and maps wave attenuation along marsh coastlines following an exponential decay.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#    This software is released under the MIT License, see LICENSE.txt.
#
### Step 1 ###########################################################
print("Step 1: Import modules and inputs parameters")
######################################################################

### 1.1 Import modules ###
print("Step 1.1: Import modules")
# import os module
import os, sys
import numpy as np
import pandas as pd
import csv
import warnings
import geopandas as gpd
import glob
# import pathlib

# import spatial analysis 
from osgeo import gdal, ogr,osr
from osgeo.gdalconst import *
gdal.AllRegister()                                      # Register all of drivers

############################### Need QGIS Module #####################
from qgis.core import *
from qgis.analysis import *
from qgis.gui import *

QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
QgsApplication.setPrefixPath(r'C:\OSGeo4W\apps\qgis-ltr', True)
qgs = QgsApplication([], False)
qgs.initQgis()
sys.path.append('C:\\OSGeo4W\\apps\\qgis\\python\\plugins')

import processing
from processing.core.Processing import Processing
Processing.initialize()
######################################################################

# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

# The target Working directory
Workspace = "Z:/CCR_data/ACTIVE/Jin_Projects/WATTE/EMBIR/Test_3/"
#Workspace = pathlib.Path.cwd().parent

# Change the current working directory
os.chdir(Workspace)

### 1.2 Input parameters and file ###
print("Step 1.2: Input parameters and file")

# Dominant wave direction North True
Wave_direction = float(150)  # North is 0 [degree] and clockwise
Wave_direction_range = 179    # Criterion for consider transects
#Wave_direction = float(input("Dominant wave direction in North True: [degree] between 0-360 type here->"))
assert 0<=Wave_direction and Wave_direction<360,"Please input the values between 0-360 degree ...{{ (>_<) }}\n"

# Input raster data
Inputspace = os.path.join(Workspace, "Input_data")
Raster_file = os.path.join(Inputspace,"Grand_Bay_2050_IM.tif")

# Input baseline delineation method
##################################################################################################################
# type 1 Method 1 Land (1): Water (0) isopleth with a moving average method
# type 2 Method 2 1-way scan based on significant wave direction
# type 3 Method 3 Baseline delineation based on a manual polyline also provide the path
# Manual_line= "Z:/CCR_data/ACTIVE/Jin_Projects/WATTE/EMBIR/Test_3/Input_data/Polyline.shp"
# Here polyline is estimated to draw from left side to right side, which affect offset direction

Baseline_delineation_method = 3
Manual_line = os.path.join(Inputspace,"Polyline.shp")
assert 1 <= Baseline_delineation_method and Baseline_delineation_method <= 3,"Please input the values between 1-3 ...{{ (>_<) }}\n"
##################################################################################################################

#Classification(s), Water, String
Input_Water_Class = "40" #Don't use 0 (no data)

#Classification(s), Other, String
Input_Other_Class= "55"

#Classification(s), Marsh, String
Input_Marsh_Class = "16,23,32"
Marsh_Class = [int(x) for x in Input_Marsh_Class.split(',')]

#Decay constant for each classification of water and marshes, String
Input_Decay_Constant = "0.021,0.030,0.090"
Decay_Constant_M = [float(x) for x in Input_Decay_Constant.split(',')]

#Transect length [m]  
Input_Transect_length = 1000

#Distance between transects [m]  
Input_Distance_Transects = 10

# Filter the wave attenuation values greater than a threshold
Input_Wave_Attenuation_Threshold = 0.02
######################################################################

#Spacing of points along transect [m]
Input_Spacing = 10                   # dx = (x_i+1−x_i)

#Interpolation_method = "IDW"          # currently IDW only

# Make process folder
Process_folder = os.path.join(Workspace, 'Process')
try:
    os.mkdir(Process_folder)
except:
    pass
# Make process folder
Outputspace = os.path.join(Workspace, 'Output_data')
try:
    os.mkdir(Outputspace)
except:
    pass
########################### Messages #################################
print("Dominant wave direction: ", Wave_direction, "[degree]")
print("Input raster dataset is " + Raster_file)
print('Water classification is ', end = " ")
print(*Input_Water_Class, sep = " ")
Input_Water_Class=int(Input_Water_Class)
print('Other classification is ', end = " ")
print(*Input_Other_Class, sep = " ")
Input_Other_Class=int(Input_Other_Class)
print('Marsh classification is ', end = " ")
print(*Marsh_Class)
print ('Decay constant for each classification of marsh is', end = " ")
print (*Decay_Constant_M)

#print(Decay_Constant_M)
print ('Transect length is', Input_Transect_length, "meter")
print ('Distance between transects is', Input_Distance_Transects, "meter")

Input_list = [Input_Water_Class] + [Input_Other_Class] + Marsh_Class
print('Input_list', Input_list)

# print ('Wave decay stopping point is', Input_Wave_decay)
print ('Spacing of points along transect is', Input_Spacing, "meter")

### Step 2 ###########################################################
print ('Step 2: Outputs')
######################################################################
Outdir=Process_folder
print ("Output Directory is ",Outdir)

### 2.1 Read raster(input) file ###
print("Step 2.1: Read raster file")

######################################################################
# Read raster(input) file
######################################################################

Rasterdata = gdal.Open(Raster_file, GA_ReadOnly)
if Rasterdata is None:
    print("Could not open " + Rasterdata)
    sys.exit(1)
print("Reading raster file (georeference, etc)")

# Coordinate system
prj = Rasterdata.GetProjection()                        # Read projection
print("Projection:", prj)

# Get raster size and band
rows = Rasterdata.RasterYSize                           # number of rows
cols = Rasterdata.RasterXSize                           # number of columns
bandnum = Rasterdata.RasterCount                        # band number
print("rows=", rows, "cols=", cols)
#print("band=", bandnum)

# Get georeference info
transform = Rasterdata.GetGeoTransform()
xOrigin = transform[0]                                  # Upperleft x
yOrigin = transform[3]                                  # Upperleft y
pixelWidth = transform[1]                               # cell size x
pixelHeight = transform[5]                              # cell size y (value is negative)
print("xOrigin=", xOrigin, "m", "yOrigin=", yOrigin, "m")
print("pixelWidth=", pixelWidth, "m", "pixelHeight=", -pixelHeight, "m") # pixelHeight is always negative

# Read the raster band
band = Rasterdata.GetRasterBand(1)
# Data type of the values
print('data type is',gdal.GetDataTypeName(band.DataType))  # Each raster file has a different data type
# Get band value info
RV = Rasterdata.GetRasterBand(1).ReadAsArray()          # raster values in the band

RV_1D = RV.reshape(-1)                                  # 1D array is needed to use set function
Input_list_data = set(RV_1D)
nodata_list = set(Input_list_data) ^ set(Input_list)    # find not matched data
nodata = list(nodata_list)

#######################################################################################################################
print('non-Input value list', nodata)                   # Need to check this
#######################################################################################################################

# Fetch metadata for the band
band.GetMetadata()

# Print only selected metadata:
print ("[ NO DATA VALUE ] = ", band.GetNoDataValue()) # check nodata
print ("[ MIN ] = ", band.GetMinimum())
print ("[ MAX ] = ", band.GetMaximum())

if band.GetMinimum() is None or band.GetMaximum() is None:
    band.ComputeStatistics(0)
    print("[ MIN ] = ", band.GetMinimum())
    print("[ MAX ] = ", band.GetMaximum())

# 2.2 Make Polygons within the domain
print("2.2 Make Polygons within the domain")

#Convert raster to polygon
Draft_Polygon = os.path.join(Outputspace,'Draft_Polygon.shp')
drv = ogr.GetDriverByName("ESRI Shapefile")
out_ds = drv.CreateDataSource(Draft_Polygon)
prj2 = osr.SpatialReference()
prj2.ImportFromWkt(Rasterdata.GetProjection())
out_layer = out_ds.CreateLayer(Draft_Polygon, srs=prj2)
New_field = ogr.FieldDefn("Value",ogr.OFTInteger)                # add Value column
out_layer.CreateField(New_field)
gdal.FPolygonize(band, None, out_layer, 0, callback=None)
out_ds.SyncToDisk()

del out_layer
del out_ds

# Classify the polygon
Domain_Polygon = os.path.join(Outputspace,'Domain_Polygon.shp')               # Polygon within Domain
Output_Land_Polygon = os.path.join(Outputspace,'Polygon_Land.shp')            # Land Polygon
Output_Water_Polygon = os.path.join(Outputspace,'Polygon_Water.shp')          # Water Polygon

# Sort the GeoDataFrame in descending order based on the Length column
gdf = gpd.read_file(Draft_Polygon)

gdf_domain = gdf[gdf['Value'] != 255]
gdf_water = gdf[gdf['Value'] == Input_Water_Class]
gdf_land = gdf[gdf['Value'] != Input_Water_Class]

# Create a new GeoDataFrame with only the selected row
gdf_domain = gpd.GeoDataFrame(gdf_domain, geometry='geometry', crs=gdf.crs)
gdf_land = gpd.GeoDataFrame(gdf_land, geometry='geometry', crs=gdf.crs)
gdf_water = gpd.GeoDataFrame(gdf_water, geometry='geometry', crs=gdf.crs)
gdf_domain.to_file(Domain_Polygon)
gdf_land.to_file(Output_Land_Polygon)
gdf_water.to_file(Output_Water_Polygon)

# Convert raster data type to Float32
NoData_value = -99999.
Raster_float = os.path.join(Workspace,'Output_data/Rasterdtype2float.tif')
gtiff_driver = gdal.GetDriverByName('GTiff')                       # Use GeoTIFF driver
out_ds = gtiff_driver.Create(Raster_float,                         # Create a output file
band.XSize,band.YSize, bandnum, gdal.GDT_Float32)
out_ds.SetProjection(prj)
out_ds.SetGeoTransform(transform)

dst_band = out_ds.GetRasterBand(1)
dst_band.WriteArray(RV)
dst_band = out_ds.GetRasterBand(1).SetNoDataValue(NoData_value)    # Exclude nodata value
dst_band = out_ds.GetRasterBand(1).ComputeStatistics(0)            # Calculate statistics for Raster pyramids (Pyramids can speed up the display of raster data)
print('maximum value in this domain is', np.max(RV))

del out_ds

# Open Raster float data and Output land and water classification raster map
Rasterfloat = gdal.Open(Raster_float, GA_ReadOnly)
if Raster_float is None:
    print('Could not open ' + Rasterfloat)
    sys.exit(1)
print("Reading raster file (georeference, etc)")

# Get band value info
AA = Rasterfloat.GetRasterBand(1).ReadAsArray()                  # raster values in the band

WL = AA
WL[WL == 0.0] = NoData_value
WL[(WL == float(Input_Water_Class))] = 0.0
WL[(WL != 0.0) & (WL != NoData_value)] = 1.0

# Output land and water classification raster map
WL_class = os.path.join(Workspace,'Output_data/WL_class.tif')
gtiff_driver = gdal.GetDriverByName('GTiff')                     # Use GeoTIFF driver
out_ds = gtiff_driver.Create(WL_class,                           # Create a output file
band.XSize,band.YSize, bandnum, gdal.GDT_Float32)
out_ds.SetProjection(prj)
out_ds.SetGeoTransform(transform)

dst_band = out_ds.GetRasterBand(1)
dst_band.WriteArray(WL)
dst_band = out_ds.GetRasterBand(1).SetNoDataValue(NoData_value)  # Exclude nodata value
dst_band = out_ds.GetRasterBand(1).ComputeStatistics(0)          # Calculate statistics for Raster pyramids (Pyramids can speed up the display of raster data)
print('maximum value of water domain probability is', np.max(WL))
del out_ds

######################################################################
# Determine the quadrant of wave direction
######################################################################

# Determine the quadrant of dominant wave direction
theta = Wave_direction * np.pi / 180.  # convert degree to radian
# print(theta)
# Set direction regime (qauadrants)
if theta > np.pi / 4. and theta <= 3 * np.pi / 4.:
    q = 41  # Quadrants 4 & 1
elif theta > 3 * np.pi / 4. and theta <= 5 * np.pi / 4.:
    q = 34  # Quadrants 3 & 4
elif theta > 5 * np.pi / 4. and theta <= 7 * np.pi / 4.:
    q = 23  # Quadrants 2 & 3
else:
    q = 12  # Quadrants 1 & 2
print("quadrants is ", q)

### Step 3 ###########################################################
print("Step 3: Make a baseline")
######################################################################

##################################################################################################################
### 3.1 Baseline delineation #######

# Method 1 Land (1): Water (0) isopleth with a moving average method
# Method 2 1-way scan based on significant wave direction
# Method 3 Baseline delineation based on a manual polyline
##################################################################################################################
if Baseline_delineation_method == 1:
    ##################################################################################################################
    # Method1 Baseline delineation based on isopleth with a moving average method
    ##################################################################################################################
    print("Baseline delineation based on isopleth with a moving average method")

    # Conduct moving average resampling (this method requires large memory. The moving windows shouldn't use a large value on your stand-alone PC)
    Moving_average_size = 31                                                # User freely change the window size
    ResampleRasterMA = os.path.join(Workspace,'Output_data/ResampleMA.tif') # Resample file name
    Contour = os.path.join(Workspace,'Output_data/Contour.shp')

    #define number of contours and range
    conNum = 1
    conList = [0.25]                                                        # Extract a contour (water/land isopleth)=[0.25,0.5, 0.75 etc]

    def make_slices(data, win_size):
    # """Return a list of slices given a window size.
    # data - two-dimensional array to get slices from
    # win_size - tuple of (rows, columns) for the moving window
    # """
        row_num = data.shape[0] - win_size[0] + 1
        col_num = data.shape[1] - win_size[1] + 1
        slices = []
        for i in range(win_size[0]):
            for j in range(win_size[1]):
                slices.append(data[i:row_num+i, j:col_num+j])
        return slices

    # Get band value info
    WLclass = gdal.Open(WL_class, GA_ReadOnly)
    indata_pre = WLclass.GetRasterBand(1).ReadAsArray()                  # raster values in the band
    indata = np.where((indata_pre == NoData_value),np.nan,indata_pre)
    slices = make_slices(indata,(Moving_average_size,Moving_average_size))
    stacked = np.dstack(slices)
    outdata = np.full(indata.shape, np.nan)
    outdata = np.zeros(indata.shape, np.float32)
    range1 = int((Moving_average_size-1)/2)                              # don't use a name "range" because range will reuse in for loops

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        outdata[range1:-range1,range1:-range1] = np.nanmean(stacked, 2)

    outdata[0:range1, :] = NoData_value
    outdata[-range1:rows, :] = NoData_value
    outdata[:, 0:range1] = NoData_value
    outdata[:, -range1:cols] = NoData_value

    gtiff_driver = gdal.GetDriverByName('GTiff')                       # Use GeoTIFF driver
    out_ds = gtiff_driver.Create(ResampleRasterMA,band.XSize,band.YSize, bandnum, gdal.GDT_Float32) # Create a output file
    out_ds.SetProjection(prj)
    out_ds.SetGeoTransform(transform)
    outdata_mod = np.where(np.isnan(outdata), NoData_value, outdata)                                # Change np.nun to numeric value
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(outdata_mod)
    out_band = out_ds.GetRasterBand(1).SetNoDataValue(NoData_value)  # Exclude nodata value
    out_band = out_ds.GetRasterBand(1).ComputeStatistics(0)          # Calculate statistics for Raster pyramids (Pyramids can speed up the display of raster data)
    print('maximum value of water domain probability is', np.max(outdata))

    del out_ds                                                       # Close the files out_ds.Destroy() or = none are also available

    # Extract contour (Land (1): Water (0) isopleth)
    ResampleRasterdata = gdal.Open(ResampleRasterMA, GA_ReadOnly)
    band =ResampleRasterdata.GetRasterBand(1)

    drv = ogr.GetDriverByName("ESRI Shapefile")                      # Set up the shapefile driver
    out_ds = drv.CreateDataSource(Contour)                           # Create a data source

    # Create a layer
    prj3 = osr.SpatialReference()
    prj3.ImportFromWkt(ResampleRasterdata.GetProjection())
    out_layer = out_ds.CreateLayer('Contour',prj3, ogr.wkbLineString)

    # Define fields of id and elev
    fieldDef = ogr.FieldDefn("ID", ogr.OFTInteger)
    out_layer.CreateField(fieldDef)
    fieldDef = ogr.FieldDefn("elev", ogr.OFTReal)
    out_layer.CreateField(fieldDef)

    # Write shapefile
    # ContourGenerate(Band srcBand, double contourInterval, double contourBase, int fixedLevelCount, int useNoData, double noDataValue,
    #                Layer dstLayer, int idField, int elevField
    gdal.ContourGenerate(band, 0, 0, conList, 1, NoData_value, out_layer, 0, 1)

    del out_ds

    # Add the length in the contour lines (you can also use shapely)
    dst_ds = ogr.Open(Contour, 1) # 0 means read-only. 1 means writeable.
    layer = dst_ds.GetLayer()
    featureCount = layer.GetFeatureCount()
    New_field = ogr.FieldDefn("Length", ogr.OFTReal)
    layer.CreateField(New_field)

    for i,feat in enumerate(layer):
        Line_length=feat.geometry().Length()
        feat.SetField("Length",Line_length)
        layer.SetFeature(feat)
    # print(feat.geometry().Length())
    del layer
    del dst_ds # Close the Shapefile

    # Sort the GeoDataFrame in descending order based on the Length column
    gdf = gpd.read_file(Contour)
    gdf = gdf.sort_values('Length', ascending=False)
    max_length_gdf = gdf.iloc[[0]]                   # Select the first row (the maximum Length value)

    # Create a new GeoDataFrame with only the selected row
    max_length_gdf = gpd.GeoDataFrame(max_length_gdf, geometry='geometry', crs=gdf.crs)
    Max_length_contour = os.path.join(Workspace, 'Output_data/Max_length_contour.shp')
    max_length_gdf.to_file(Max_length_contour)

    #Create contour line points
    Output_point=os.path.join(Outputspace,'Contour_points.shp')
    processing.run("native:pointsalonglines", {'INPUT':Max_length_contour,'DISTANCE':Input_Distance_Transects,'START_OFFSET':0,'END_OFFSET':0,'OUTPUT':Output_point})

    Output_point_csv=os.path.join(Process_folder,'Contour_Points.csv')
    gdal.VectorTranslate(Output_point_csv, Output_point, format='CSV', layerCreationOptions=['GEOMETRY=AS_WKT'])

    with open(Output_point_csv, "r", newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)  # read the header row
        contents = list(csvreader)    # read the remaining rows

    dst_ds = ogr.Open(Output_point, 0)  # 0 means read-only. 1 means writeable.
    if dst_ds is None:
        sys.exit('Could not open {0}.'.format(dst_ds))
    layer = dst_ds.GetLayer()
    featureCount = layer.GetFeatureCount()
    fields = ["East", "North"]                                              # write field names
    header.extend(fields)
    for i, feat in enumerate(layer):
    # feat = layer.GetFeature(i)
        wkt_point = feat.GetGeometryRef().ExportToWkt()
        point = ogr.CreateGeometryFromWkt(wkt_point)
        x, y = point.GetX(), point.GetY()
        contents[i].append(x)
        contents[i].append(y)

    contents.reverse()                                                      # reorder the contents When polyline direction is East to West side
    with open(Output_point_csv, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(header)
        csvwriter.writerows(contents)                                       # writing data

    del dst_ds                                                              # Close file

    df = pd.read_csv(Output_point_csv)                                      # Using pandas in this example
    m1=df.shape[0]
    moving_value=5
    transiton_points=3
    file="Contour_Points_MA"+str(int(moving_value))
    file_name=file
    #file = df.loc[:,['East', 'North']].rolling(moving_value,center=True).mean()                ## NaN values will be replaced by original values
    file = df.loc[:,['East', 'North']].rolling(moving_value,min_periods=moving_value-transiton_points*2,center=True).mean()                ## NaN values will be replaced by original values
    file.loc[:int(np.floor(moving_value/2)-transiton_points),['East', 'North']]=df.loc[:int(np.floor(moving_value/2-transiton_points)),['East', 'North']]
    file.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]=df.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]

    save_name=file_name+'.csv'
    save_name2=os.path.join(Process_folder,save_name)
    file.to_csv(save_name2, index=False)

    # Point to polyline
    line = ogr.Geometry(ogr.wkbLineString)
    for i in df.index:
        line.AddPoint(float(file.East[i]), float(file.North[i]))

    Smoothed_polyline = os.path.join(Outputspace,'Smoothed_line_Contour.shp')       # output file name and location
    drv = ogr.GetDriverByName("ESRI Shapefile")                      # Set up the shapefile driver
    out_ds = drv.CreateDataSource(Smoothed_polyline)                 # Create a data source
    prj5 = osr.SpatialReference()
    prj5.ImportFromWkt(Rasterdata.GetProjection())
    # Create a layer
    out_layer = out_ds.CreateLayer('Smoothed_polyline', prj5, ogr.wkbLineString)  # Create a layer

    # Add an ID field in the layer
    idField = ogr.FieldDefn('id', ogr.OFTInteger)                    # Add a column 'id'
    out_layer.CreateField(idField)

    # Create a feature and set geometry
    feature_Defn = out_layer.GetLayerDefn()                          # Create a blank(dummy) feature
    feature = ogr.Feature(feature_Defn)
    feature.SetGeometry(line)                                        # Set polygon
    feature.SetField('id', 1)                                        # Add id number
    out_layer.CreateFeature(feature)                                 # Save the feature in the layer

    del feature                                                      # destroy/delete feature
    del out_ds                                                       # Close file

    ##################################################################################################################
    # Method2 Baseline delineation based on a 1-way scan based on significant wave direction
    ##################################################################################################################
elif Baseline_delineation_method == 2:
    print("Baseline delineation based on a 1-way scan based on significant wave direction")

    ### Find interface (shoreline)points ###
    print("Find interface (shoreline) points")

    #Determine scanning direction
    if q == 34:
        RV_ad = np.flip(RV, 0)  # bottom(South) to top(North) Quadrants 3 & 4
    elif q == 12:
        RV_ad = RV              # top(North) to bottom(South) Quadrants 1 & 2
    elif q == 41:
        RV_ad = np.flip(RV, 1)  # right(East) to left(West) Quadrants 4 & 1
    else:
        RV_ad = RV              # left(West) to right(East) Quadrants 2 & 3

    ### Output band value (raw) ###
    print("Output raster values in csv file")

    csv_out = os.path.join(Process_folder,'raster_values_raw.csv')
    with open(csv_out, "w+", newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerows(RV)
        csvfile.close()

    # Find interface between water and land)
    print("Find interface points between water and land")

    if q == 12 or q == 34:  # Quadrants 1 & 2 or # Quadrants 3 & 4
        fp = np.zeros((cols, 2))  # first point (column, row)
        j = 0
        for j in range(cols):
            #    print(j)
            i = 0
            for i in range(rows):
                if RV_ad[i, j] == 0 or RV_ad[i, j] == nodata or RV_ad[i, j] == Input_Water_Class:
                    i += 1
                else:
                    # print(i+1,j)
                    break
            fp[j, 0] = j
            fp[j, 1] = i
            j += 1
    # print (fp)
    else: # Quadrants 2 & 3 or # Quadrants 4 & 1
        fp = np.zeros((rows, 2))  # first point (column, row)
        i = 0
        for i in range(rows):
            #    print(j)
            j = 0
            for j in range(cols):
                if RV_ad[i, j] == 0 or RV_ad[i, j] == nodata or RV_ad[i, j] == Input_Water_Class:
                    j += 1
                else:
                    # print(i+1,j)
                    break
            fp[i, 0] = i
            fp[i, 1] = j
            i += 1
    # print (fp)

    csv_out = os.path.join(Process_folder,'fp.csv')
    with open(csv_out, "w+", newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerows(fp)                                                  # writing data
        csvfile.close()

    ### Compute coordinates of the points ###
    print("Compute and save coordinates of the points (x,y)")

    if q == 12:
        xy = np.zeros((cols, 2))
        for j in range(cols):

            xy[j, 0] = xOrigin + pixelWidth / 2 + fp[j, 0] * pixelWidth
            xy[j, 1] = yOrigin - (fp[j, 1]) * (-pixelHeight)         # pixel height is negative

    elif q == 34:
        xy = np.zeros((cols, 2))
        for j in range(cols):

            xy[j, 0] = xOrigin + pixelWidth/2 + fp[j, 0] * pixelWidth
            xy[j, 1] = yOrigin - (rows - fp[j, 1]) * (-pixelHeight)  # pixel height is negative

    elif q == 23:
        xy = np.zeros((rows, 2))
        for i in range(rows):

            xy[i, 0] = xOrigin + fp[i, 1] * pixelWidth
            xy[i, 1] = yOrigin -(-pixelHeight)/2 - (fp[i, 0]) * (-pixelHeight)  # pixel height is negative

    else:
        xy = np.zeros((rows, 2))
        for i in range(rows):

            xy[i, 0] = xOrigin + (cols-fp[i, 1]) * pixelWidth
            xy[i, 1] = yOrigin -(-pixelHeight)/2 - (fp[i, 0]) * (-pixelHeight)  # pixel height is negative

    csv_out = os.path.join(Process_folder,'coordinate_W.csv')
    with open(csv_out, "w+", newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        fields = ["East", "North"]                                              # write field names
        csvwriter.writerow(fields)
        csvwriter.writerows(xy)                                                 # writing data
        csvfile.close()

    Input_point_W = os.path.join(Outputspace,'Input_points_W.shp')
    drv = ogr.GetDriverByName("ESRI Shapefile")
    out_ds = drv.CreateDataSource(Input_point_W)
    prj4 = osr.SpatialReference()
    prj4.ImportFromWkt(Rasterdata.GetProjection())
    out_layer = out_ds.CreateLayer('Input_Points', prj4, ogr.wkbPoint)
    layer_defn = out_layer.GetLayerDefn()                                       # gets parameters of the current shapefile
    index = 0

    # reading the points from csv file
    with open(csv_out, 'r') as csvfile:
        readerDict = csv.DictReader(csvfile, delimiter=',')  # this code uses the tab as delimiter
        for field in readerDict.fieldnames:
            new_field = ogr.FieldDefn(field, ogr.OFTString)  # Create a new field with the content of header
            out_layer.CreateField(new_field)
        for row in readerDict:
            #print(row['East'], row['North'])
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(float(row['East']), float(row['North']))  # We have 'East' and 'North' values as Strings, so convert into float
            feature = ogr.Feature(layer_defn)
            feature.SetGeometry(point)                               # set the coordinates
            feature.SetFID(index)
            for field in readerDict.fieldnames:
                i = feature.GetFieldIndex(field)
                feature.SetField(i, row[field])
            out_layer.CreateFeature(feature)
            index += 1

    del out_ds                                                       # Save and close the data source

    dst_ds = ogr.Open(Input_point_W, 1) # 0 means read-only. 1 means writeable.
    layer = dst_ds.GetLayer()
    featureCount = layer.GetFeatureCount()

    New_field = ogr.FieldDefn("FID", ogr.OFTInteger)
    layer.CreateField(New_field)

    for i,feat in enumerate(layer):
        feat.SetField('FID', i)
        layer.SetFeature(feat)

    del layer
    del dst_ds # Close the Shapefile

    ############################### Potential future modification #############################
    # #Extract Points within water region
    Output_point = os.path.join(Outputspace,'Points.shp')
    processing.run("native:extractbylocation", {'INPUT':Input_point_W,'PREDICATE':[0],'INTERSECT':Output_Water_Polygon,'OUTPUT':Output_point})

    Output_point_csv = os.path.join(Process_folder,'Shoreline_Points_W.csv')
    gdal.VectorTranslate(Output_point_csv, Output_point, format='CSV', layerCreationOptions=['GEOMETRY=AS_WKT'])

    df = pd.read_csv(Output_point_csv)                                           # Using pandas in this example
    m1 = df.shape[0]
    moving_value = 5
    transiton_points = 3
    file = "Shoreline_Points_W_MA"+str(int(moving_value))
    file_name = file
    #file = df.loc[:,['East', 'North']].rolling(moving_value,center=True).mean()                ## NaN values will be replaced by original values
    file = df.loc[:,['East', 'North']].rolling(moving_value,min_periods=moving_value-transiton_points*2,center=True).mean()                ## NaN values will be replaced by original values
    file.loc[:int(np.floor(moving_value/2)-transiton_points),['East', 'North']]=df.loc[:int(np.floor(moving_value/2-transiton_points)),['East', 'North']]
    file.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]=df.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]

    save_name = file_name+'.csv'
    save_name2 = os.path.join(Process_folder,save_name)
    file.to_csv(save_name2, index=False)

    # Point to polyline
    line = ogr.Geometry(ogr.wkbLineString)
    for i in df.index:
        line.AddPoint(float(file.East[i]), float(file.North[i]))

    Smoothed_polyline = os.path.join(Outputspace,'Smoothed_line_W.shp')       # output file name and location
    drv = ogr.GetDriverByName("ESRI Shapefile")                               # Set up the shapefile driver
    out_ds = drv.CreateDataSource(Smoothed_polyline)                          # Create a data source
    prj5 = osr.SpatialReference()
    prj5.ImportFromWkt(Rasterdata.GetProjection())
    # Create a layer
    out_layer = out_ds.CreateLayer('Smoothed_polyline', prj5, ogr.wkbLineString)  # Create a layer

    # Add an ID field in the layer
    idField = ogr.FieldDefn('id', ogr.OFTInteger)                    # Add a column 'id'
    out_layer.CreateField(idField)

    # Create a feature and set geometry
    feature_Defn = out_layer.GetLayerDefn()                          # Create a blank(dummy) feature
    feature = ogr.Feature(feature_Defn)
    feature.SetGeometry(line)                                        # Set polygon
    feature.SetField('id', 1)                                        # Add id number
    out_layer.CreateFeature(feature)                                 # Save the feature in the layer

    del feature                                                      # destroy/delete feature
    del out_ds                                                       # Close file

    ##################################################################################################################
    # Method3 Baseline delineation based on a manual input
    ##################################################################################################################
else:
    print("Baseline delineation based on a manual polyline")
    #Create manual line points
    Output_point=os.path.join(Outputspace,'Manual_points.shp')
    processing.run("native:pointsalonglines", {'INPUT':Manual_line,'DISTANCE':Input_Distance_Transects,'START_OFFSET':0,'END_OFFSET':0,'OUTPUT':Output_point})

    Output_point_csv=os.path.join(Process_folder,'Manual_line_Points.csv')
    gdal.VectorTranslate(Output_point_csv, Output_point, format='CSV', layerCreationOptions=['GEOMETRY=AS_WKT'])

    with open(Output_point_csv, "r", newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)  # read the header row
        contents = list(csvreader)    # read the remaining rows

    dst_ds = ogr.Open(Output_point, 0)  # 0 means read-only. 1 means writeable.
    if dst_ds is None:
        sys.exit('Could not open {0}.'.format(dst_ds))
    layer = dst_ds.GetLayer()
    featureCount = layer.GetFeatureCount()
    fields = ["East", "North"]                                              # write field names
    header.extend(fields)
    for i, feat in enumerate(layer):
    # feat = layer.GetFeature(i)
        wkt_point = feat.GetGeometryRef().ExportToWkt()
        point = ogr.CreateGeometryFromWkt(wkt_point)
        x, y = point.GetX(), point.GetY()
        contents[i].append(x)
        contents[i].append(y)

    #contents.reverse()                                                     # reorder the contents When polyline direction is East to West side
    with open(Output_point_csv, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(header)
        csvwriter.writerows(contents)                                       # writing data

    del dst_ds                                                              # Close file

    df = pd.read_csv(Output_point_csv)                                      # Using pandas in this example
    m1 = df.shape[0]
    moving_value = 5
    transiton_points = 3
    file = "Manual_line_Points_MA"+str(int(moving_value))
    file_name = file
    #file = df.loc[:,['East', 'North']].rolling(moving_value,center=True).mean()                ## NaN values will be replaced by original values
    file = df.loc[:,['East', 'North']].rolling(moving_value,min_periods=moving_value-transiton_points*2,center=True).mean()                ## NaN values will be replaced by original values
    file.loc[:int(np.floor(moving_value/2)-transiton_points),['East', 'North']]=df.loc[:int(np.floor(moving_value/2-transiton_points)),['East', 'North']]
    file.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]=df.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]

    save_name = file_name+'.csv'
    save_name2 = os.path.join(Process_folder,save_name)
    file.to_csv(save_name2, index=False)

    #Point to polyline
    line = ogr.Geometry(ogr.wkbLineString)
    for i in df.index:
        line.AddPoint(float(file.East[i]), float(file.North[i]))

    Smoothed_polyline = os.path.join(Outputspace,'Smoothed_line_Manual.shp')   # output file name and location
    drv = ogr.GetDriverByName("ESRI Shapefile")                                # Set up the shapefile driver
    out_ds = drv.CreateDataSource(Smoothed_polyline)                           # Create a data source
    prj5 = osr.SpatialReference()
    prj5.ImportFromWkt(Rasterdata.GetProjection())
    # Create a layer
    out_layer = out_ds.CreateLayer('Smoothed_polyline', prj5, ogr.wkbLineString)  # Create a layer

    # Add an ID field in the layer
    idField = ogr.FieldDefn('id', ogr.OFTInteger)                    # Add a column 'id'
    out_layer.CreateField(idField)

    # Create a feature and set geometry
    feature_Defn = out_layer.GetLayerDefn()                          # Create a blank(dummy) feature
    feature = ogr.Feature(feature_Defn)
    feature.SetGeometry(line)                                        # Set polygon
    feature.SetField('id', 1)                                        # Add id number
    out_layer.CreateFeature(feature)                                 # Save the feature in the layer

    del feature                                                      # destroy/delete feature
    del out_ds                                                       # Close file

###########################################################################################
# Offset to offshore direction and Make a baseline
############################### Potential future modification #############################

Offset_line_pre = os.path.join(Outputspace,'Offset_Line_pre.shp')
#Offset_line_old = os.path.join(Outputspace,'Offset_Line_qgis.shp')
Dist = (30*((abs(pixelWidth)+abs(pixelHeight))/2))  # 30* average cell size
if q == 23 or q == 34:
    Dist_offset = -Dist
else:
    Dist_offset = Dist
processing.run("gdal:offsetcurve", {'INPUT':Smoothed_polyline,'GEOMETRY':'geometry','DISTANCE':Dist_offset,'OPTIONS':'','OUTPUT':Offset_line_pre}) #this one is better
print('done offset')

# Read the offset shapefile
gdf = gpd.read_file(Offset_line_pre)
# Explode the geometry column to separate the disconnected lines
gdf_exploded = gdf.explode(column='geometry', index_parts=True)

# Save the exploded data to a new shapefile
Offset_line = os.path.join(Outputspace,'Offset_Line.shp')
gdf_exploded.to_file(Offset_line)

# Add the length in the offset line (you can also use shapely)
dst_ds = ogr.Open(Offset_line, 1)  # 0 means read-only. 1 means writeable.
layer = dst_ds.GetLayer()
featureCount = layer.GetFeatureCount()
New_field = ogr.FieldDefn("Length", ogr.OFTReal)
layer.CreateField(New_field)

for i, feat in enumerate(layer):
    Line_length = feat.geometry().Length()
    feat.SetField("Length", Line_length)
    layer.SetFeature(feat)
# print(feat.geometry().Length())
del layer
del dst_ds  # Close the Shapefile

# Sort the GeoDataFrame in descending order based on the Length column
gdf = gpd.read_file(Offset_line)
gdf = gdf.sort_values('Length', ascending=False)
max_length_gdf = gdf.iloc[[0]]  # Select the first row (the maximum Length value)

# Create a new GeoDataFrame with only the selected row
max_length_gdf = gpd.GeoDataFrame(max_length_gdf, geometry='geometry', crs=gdf.crs)
Baseline_pre = os.path.join(Workspace, 'Output_data/Baseline_pre.shp')
max_length_gdf.to_file(Baseline_pre)
del max_length_gdf                                                          # geopandas is still open or blocked to delete

# Create baseline_pre points
Baseline_pre_points=os.path.join(Outputspace,'Baseline_pre_points.shp')
processing.run("native:pointsalonglines", {'INPUT':Baseline_pre,'DISTANCE':Input_Distance_Transects,'START_OFFSET':0,'END_OFFSET':0,'OUTPUT':Baseline_pre_points})
# Note: the point oder is East (start point) to West (end point) now due to offset command.

Output_point_csv=os.path.join(Process_folder,'Baseline_pre_Points.csv')
gdal.VectorTranslate(Output_point_csv, Baseline_pre_points, format='CSV', layerCreationOptions=['GEOMETRY=AS_WKT'])

with open(Output_point_csv, "r", newline='') as csvfile:
    csvreader = csv.reader(csvfile)
    header = next(csvreader)  # read the header row
    contents = list(csvreader)    # read the remaining rows

dst_ds = ogr.Open(Baseline_pre_points, 0)  # 0 means read-only. 1 means writeable.
if dst_ds is None:
    sys.exit('Could not open {0}.'.format(dst_ds))
layer = dst_ds.GetLayer()
featureCount = layer.GetFeatureCount()
fields = ["East", "North"]                                              # write field names
header.extend(fields)
for i, feat in enumerate(layer):
# feat = layer.GetFeature(i)
    wkt_point = feat.GetGeometryRef().ExportToWkt()
    point = ogr.CreateGeometryFromWkt(wkt_point)
    x, y = point.GetX(), point.GetY()
    contents[i].append(x)
    contents[i].append(y)

#contents.reverse()                                                     # reorder the contents When polyline direction is East to West side
with open(Output_point_csv, "w", newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header)
    csvwriter.writerows(contents)                                       # writing data

del layer
del dst_ds                                                              # Close file

df = pd.read_csv(Output_point_csv)                                      # Using pandas in this example
m1 = df.shape[0]
moving_value = 41
transiton_points = 3
file = "Baseline"
file_name = file
#file = df.loc[:,['East', 'North']].rolling(moving_value,center=True).mean()                ## NaN values will be replaced by original values
file = df.loc[:,['East', 'North']].rolling(moving_value,min_periods=moving_value-transiton_points*2,center=True).mean()                ## NaN values will be replaced by original values
file.loc[:int(np.floor(moving_value/2)-transiton_points),['East', 'North']]=df.loc[:int(np.floor(moving_value/2-transiton_points)),['East', 'North']]
file.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]=df.loc[m1-int(np.floor(moving_value/2+transiton_points)):,['East', 'North']]

save_name = file_name+'.csv'
save_name2 = os.path.join(Process_folder,save_name)
file.to_csv(save_name2, index=False)

#Point to polyline
line = ogr.Geometry(ogr.wkbLineString)
for i in df.index:
    line.AddPoint(float(file.East[i]), float(file.North[i]))

Baseline = os.path.join(Outputspace,'Baseline.shp')   # output file name and location
drv = ogr.GetDriverByName("ESRI Shapefile")                                # Set up the shapefile driver
out_ds = drv.CreateDataSource(Baseline)                                    # Create a data source
prj5 = osr.SpatialReference()
prj5.ImportFromWkt(Rasterdata.GetProjection())
# Create a layer
out_layer = out_ds.CreateLayer('Baseline', prj5, ogr.wkbLineString)  # Create a layer

# Add an ID field in the layer
idField = ogr.FieldDefn('id', ogr.OFTInteger)                    # Add a column 'id'
out_layer.CreateField(idField)

# Create a feature and set geometry
feature_Defn = out_layer.GetLayerDefn()                          # Create a blank(dummy) feature
feature = ogr.Feature(feature_Defn)
feature.SetGeometry(line)                                        # Set polygon
feature.SetField('id', 1)                                        # Add id number
out_layer.CreateFeature(feature)                                 # Save the feature in the layer

del feature                                                      # destroy/delete feature
del out_ds                                                       # Close file
print('done baseline')

# #Create baseline points
Baseline_points=os.path.join(Outputspace,'Baseline_points.shp')
processing.run("native:pointsalonglines", {'INPUT':Baseline,'DISTANCE':Input_Distance_Transects,'START_OFFSET':0,'END_OFFSET':0,'OUTPUT':Baseline_points})
# Note: the point oder is East (start point) to West (end point) now due to offset command.

### Step 4 ###########################################################
print("Step 4: Make transect")
######################################################################
### 4.1 Compute coordinates of the mid/end points for two transects ###
print("Step 4.1: Compute coordinates of the mid/end points for two transects")

dis = Input_Transect_length+1                                         ### Avoid python computation errors, add aditional distance e.g. +1 m
NP = int(np.ceil(dis/Input_Spacing))                                  # number of points on the transects, don't use np becasue we alreaqdy use np as numpy shortcut
fn = Baseline_points                                                  # vector file name
dst_ds = ogr.Open(fn, 0)                                              # 0 means read-only. 1 means writeable.
if dst_ds is None:
    sys.exit('Could not open {0}.'.format(fn))
layer = dst_ds.GetLayer()
featureCount = layer.GetFeatureCount()

print("Number of features in %s: %d" % (os.path.basename(fn),featureCount))

fp = np.zeros((featureCount, 10))                                      # feature input (column, row)
for i,feat in enumerate(layer):
    #feat = layer.GetFeature(i)
    point = feat.GetGeometryRef()
    x, y = point.GetX(), point.GetY()
    azimuth = feat["angle"]                                            # Get the angle from the shapefile feat[5] is feat["angle"]
    az1 = azimuth + 90                                                 # azimuth1 positive direction
    if az1 > 360:
        az1 -= 360
    else:
        az1 = az1

    az2 = azimuth - 90                                                 # azimuth2 negative direction
    if az2 < 0:
        az2 += 360
    else:
        az2 = az2

    fp[i, 0] = x
    fp[i, 1] = y
    fp[i, 2] = az1
    fp[i, 3] = az2
    i += 1

mx = np.zeros((featureCount-1, 1))  # mid point x  (column, row)
my = np.zeros((featureCount-1, 1))  # mid point y (column, row)
#print((fp[i, 0] + fp[i+1, 0])/2)

for i in range(0, featureCount-1):
    mx[i,0] = (fp[i, 0] + fp[i+1, 0])/2
    my[i,0] = (fp[i, 1] + fp[i+1, 1])/2
    fp[i,0] = mx[i,0]
    fp[i,1] = my[i,0]

    # Convert to Cartesian coordinate
    if -(fp[i,2]-90) < 0:
        az1_Cartesian = 360-(fp[i,2]-90)
    else:
        az1_Cartesian = -(fp[i,2] - 90)

    fp[i,4] = az1_Cartesian

    dx2 = dis*np.cos(az1_Cartesian* np.pi / 180.)
    x2 = mx[i,0]+dx2
    dy2 = dis*np.sin(az1_Cartesian* np.pi / 180.)
    y2 = my[i,0]+dy2

    fp[i,6] = x2
    fp[i,7] = y2
    if -(fp[i,3]-90) < 0:
        az2_Cartesian = 360-(fp[i,3]-90)
    else:
        az2_Cartesian = -(fp[i,3] - 90)

    fp[i,5] = az2_Cartesian

    dx3 = dis*np.cos(az2_Cartesian* np.pi / 180.)
    x3 = mx[i,0]+dx3
    dy3 = dis*np.sin(az2_Cartesian* np.pi / 180.)
    y3 = my[i,0]+dy3

    fp[i,8] = x3
    fp[i,9] = y3

    #print(az1_Cartesian)
    i += 1

fp.resize(featureCount-1, 10)      # delete last row
#print(fp.shape)
#print(fp)

# Output coordinates of the points as a csv file
csv_out = os.path.join(Process_folder,'Azimuth.csv')
with open(csv_out, "w+", newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    fields = ["mx", "my","az1","az2","az1_c","az2_c","x2","y2","x3","y3"]     # writing field names
    csvwriter.writerow(fields)
    csvwriter.writerows(fp)                                                   # writing data
    csvfile.close()

### 4.2 Make two transects
print("Step 4.2: Make two transects")
# Read Azimuth.csv
df = pd.read_csv(csv_out)                                                     # using pandas

paths_a = ogr.Geometry(ogr.wkbMultiLineString)                                # azimuth1 on the right side facing forward
paths_b = ogr.Geometry(ogr.wkbMultiLineString)                                # azimuth2 on the left side facing forward

###### paths_a and paths_b are sensitively affected by Wave_direction_range ##########################################################################
###### Users can switch the logical statemenet depending on the situations ###########################
for i in range(df.shape[0]):
    if abs(Wave_direction-df.iloc[i,3]) <= Wave_direction_range: #abs(Wave_direction-df.iloc[i,3]) != None:
        path_a = ogr.Geometry(ogr.wkbLineString)
        path_a.AddPoint(df.iloc[i,0], df.iloc[i,1])
        path_a.AddPoint(df.iloc[i,6], df.iloc[i,7])
        paths_a.AddGeometry(path_a)

    if abs(Wave_direction-df.iloc[i,2]) <= Wave_direction_range: #abs(Wave_direction-df.iloc[i,2]) != None:
        path_b = ogr.Geometry(ogr.wkbLineString)
        path_b.AddPoint(df.iloc[i,0], df.iloc[i,1])
        path_b.AddPoint(df.iloc[i,8], df.iloc[i,9])
        paths_b.AddGeometry(path_b)
#######################################################################################################################################################
# azimuth1 on the right side facing forward
Transect_a = os.path.join(Outputspace,'Transect_a.shp')
drv = ogr.GetDriverByName("ESRI Shapefile")                         # Set up the shapefile driver
out_ds = drv.CreateDataSource(Transect_a)                           # Create the data source
prj6 = osr.SpatialReference()
prj6.ImportFromWkt(Rasterdata.GetProjection())

layer = out_ds.CreateLayer("Transect_a", prj6, ogr.wkbLineString)   # Create a layer
idField = ogr.FieldDefn("id", ogr.OFTInteger)                       # Add an ID field
layer.CreateField(idField)

for i,Ft in enumerate(paths_a):
    featureDefn = layer.GetLayerDefn()                              # Create the feature and set values
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(Ft)
    feature.SetField("id", i)
    layer.CreateFeature(feature)
    del feature
na_transect = layer.GetFeatureCount()                               # Number of transect1
del out_ds                                                          # Save and close the data source

# azimuth2 on the left side facing forward
Transect_b = os.path.join(Outputspace,'Transect_b.shp')
drv = ogr.GetDriverByName("ESRI Shapefile")                         # Set up the shapefile driver
out_ds = drv.CreateDataSource(Transect_b)                           # Create the data source

layer = out_ds.CreateLayer("Transect_b", prj6, ogr.wkbLineString)   # Create a layer
idField = ogr.FieldDefn("id", ogr.OFTInteger)                       # Add an ID field
layer.CreateField(idField)

for i,xb in enumerate(paths_b):
    # Create the feature and set values
    featureDefn = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(xb)
    feature.SetField("id", i)
    layer.CreateFeature(feature)
    del feature
nb_transect = layer.GetFeatureCount()                              # Number of transect2
del out_ds                                                         # Save and close the data source

print("Number of transect1 =",na_transect,",Number of transect2 =",nb_transect)
print("Need to check the number of transect1 or 2 is zero")
print("Successfully made transects (* ^ ω ^)")

### 4.3 Convert polyline to points and read raster values ###
print("Step 4.3 Convert polyline to points and copy raster values")

###Select correct transect ###
if na_transect >= nb_transect:
    transect = Transect_a
    n_transect = na_transect
else:
    transect = Transect_b
    n_transect = nb_transect

# Create transect points
Transect_points = os.path.join(Outputspace,'Transect_points.shp')
processing.run("native:pointsalonglines",{ 'DISTANCE' : Input_Spacing, 'END_OFFSET' : 0, 'INPUT' : transect, 'OUTPUT' : Transect_points, 'START_OFFSET' : 0 })

# Copy the raster values
RV_transect_points=os.path.join(Outputspace,'RV_transect_points.shp')
processing.run("native:rastersampling",{'COLUMN_PREFIX': 'Raster_Val', 'INPUT': Transect_points, 'OUTPUT': RV_transect_points, 'RASTERCOPY': Raster_file })

### Step 5 ###########################################################
print("Step 5: Calculate wave attenuation")
######################################################################
### 5.1 Add decay constants ###
print("Step 5.1: Add decay constants")

# Read raster values
fn2 = RV_transect_points
dst_ds = ogr.Open(fn2, 1) # 0 means read-only. 1 means writeable.
if dst_ds is None:
    sys.exit('Could not open {0}.'.format(fn2))
layer = dst_ds.GetLayer()
featureCount = layer.GetFeatureCount()
print("Number of features in %s: %d" % (os.path.basename(fn2),featureCount))

# Get field names
featureDefn = layer.GetLayerDefn()
Field_names = [featureDefn.GetFieldDefn(i).GetName() for i in range(featureDefn.GetFieldCount())]
#print(Field_names)

## Add decay constant
New_field = ogr.FieldDefn('DecayC', ogr.OFTReal)
layer.CreateField(New_field)

for feat in layer:
    class_val = feat.GetField('Raster_Val')
    if class_val in Marsh_Class:
        decay_val = Decay_Constant_M[Marsh_Class.index(class_val)]
        feat.SetField('DecayC', decay_val)
    elif class_val == Input_Water_Class:
        feat.SetField('DecayC', 0.0)
    else:
        feat.SetField('DecayC', 9999)
    layer.SetFeature(feat)

del dst_ds # Close the Shapefile

### 5.2 Calculate wave attenuation ###
print("Step 5.2 calculate wave attenuation")

# Calculate wave attenuation
dst_ds = ogr.Open(fn2, 0) # 0 means read-only. 1 means writeable.
layer = dst_ds.GetLayer()
featureCount = layer.GetFeatureCount()

aa = []                                               ## processing data for decay constant
for feat in layer:
    DC = feat.GetField('DecayC')
    aa.append(DC)
aa = np.array(aa)

# convert to size to caluculate the wave attenuation value on each transect
k_const=np.reshape(aa,(n_transect, NP))               #n_transect: number of transect, NP: number of points on the transect

csv_out=os.path.join(Process_folder,'Decay_Const.csv')
with open(csv_out, "w+", newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerows(k_const)
    csvfile.close()

WT = np.zeros((k_const.shape[0],k_const.shape[1]),'d')                 # output double

for i in range(k_const.shape[0]):
    WT_pre = 1
    for j in range(k_const.shape[1]):
        k = k_const[i,j]
        if k == 9999:                                                  # 9999 means land region and no water propagate downstream
            WT[i,j] = np.nan
            WT_pre = np.nan
        else:
            WT_Temp = WT_pre*np.exp(-k*int(float(Input_Spacing)))
            WT[i,j] = WT_Temp
            WT_pre = WT_Temp

del dst_ds # Close the Shapefile

csv_out=os.path.join(Process_folder,'WT.csv')
with open(csv_out, "w+", newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerows(WT)
    csvfile.close()

# Reshape the size
WT_reshape = np.reshape(WT,(featureCount,1))

# Add Wave attenuation layer
dst_ds = ogr.Open(fn2, 1) # 0 means read-only. 1 means writeable.
layer = dst_ds.GetLayer()
featureCount = layer.GetFeatureCount()

New_field = ogr.FieldDefn('WT', ogr.OFTReal)
layer.CreateField(New_field)

if featureCount == WT_reshape.shape[0]:
    for i,feat in enumerate(layer):
        feat.SetField('WT', float(WT_reshape[i]))
        layer.SetFeature(feat)

New_field = ogr.FieldDefn('x', ogr.OFTReal)
layer.CreateField(New_field)

New_field = ogr.FieldDefn('y', ogr.OFTReal)
layer.CreateField(New_field)

for i,feat in enumerate(layer):
    #feat = layer.GetFeature(i)
    point = feat.GetGeometryRef()
    feat.SetField('x', point.GetX())
    feat.SetField('y', point.GetY())
    layer.SetFeature(feat)

    # Set the feature geometry using the point
    feat.SetGeometry(point)
    # Destroy the feature to free resources
    feat.Destroy()

del dst_ds # Close the Shapefile

### 5.3 delete unnecessary wave attenuation points (small wt value affects the results.) ###
print("Step 5.3 delete unnecessary wave attenuation points")

### 1 Filter by Input_Wave_Attenuation_Threshold using geopandas
Input_Wave_Attenuation_Threshold

gdf = gpd.read_file(fn2)

# Create a new geodataframe with only records where WT > 0.02
gdf_new = gdf[gdf['WT'] > Input_Wave_Attenuation_Threshold]

# Save the new shapefile
fn3=os.path.join(Outputspace,'RV_transect_points_threshold.shp')
gdf_new.to_file(fn3)
del gdf_new                                                       # geopandas is still open or blocked to delete

### 2 Create a Fishnet ###

spacing = round(float(Input_Spacing)*3/4,2)

dst_ds = ogr.Open(fn3, 0)
layer = dst_ds.GetLayer()
extent = layer.GetExtent()
print("extent of points id",extent)
del dst_ds

x_width = extent[1]-extent[0]
y_height = extent[3]-extent[2]
transform3 = (extent[0]-spacing/2, spacing, 0.0, extent[3]+spacing/2, 0.0, -spacing)

row_num = int(np.ceil(y_height/spacing)+1)
col_num = int(np.ceil(x_width/spacing)+1)

AA = np.zeros((row_num,col_num), dtype=int)                      # Duplication .copy function is immutable (NOT Deep copy)
id = 0
for i in range(row_num):
    for j in range(col_num):
        id += 1
        AA[i,j] = id
AA.shape

Fishnet_file = os.path.join(Outputspace,'Fishnet.tif')           # output filename
gtiff_driver = gdal.GetDriverByName('GTiff')                     # Use GeoTIFF driver
out_ds = gtiff_driver.Create(Fishnet_file,                       # Create a output file
col_num, row_num, bandnum, GDT_UInt32)

prj7 = Rasterdata.GetProjection()
out_ds.SetProjection(prj7)
out_ds.SetGeoTransform(transform3)                               # determine the position

out_band = out_ds.GetRasterBand(1)
out_band.WriteArray(AA)
out_band = out_ds.GetRasterBand(1).ComputeStatistics(0)          # calculate statistics for Raster pyramids (Pyramids can speed up the display of raster data)
print('maximum id in this fishnet is', np.max(AA))

del out_ds
RV_transect_points_fishnet = os.path.join(Outputspace,'RV_transect_points_fishnet.shp')

# overwrite  RV_transect_points
processing.run("native:rastersampling", {'INPUT':fn3,'RASTERCOPY':Fishnet_file,'COLUMN_PREFIX':'FishnetID_','OUTPUT':RV_transect_points_fishnet})
##################################################################################################
RV_transect_points_fishnet_csv=os.path.join(Process_folder,'RV_transect_points_fishnet.csv')
gdal.VectorTranslate(RV_transect_points_fishnet_csv, RV_transect_points_fishnet, format='CSV', layerCreationOptions=['GEOMETRY=AS_WKT'])

df = pd.read_csv(RV_transect_points_fishnet_csv)                                           # Using pandas in this example
# print (df)
# print (df.dtypes)

# Filtering the WT points (can be replaced by geopandas)
df2 = df.sort_values(by=['FishnetID_','WT'],ascending=False).groupby('FishnetID_').head(1) # only select maximum value
df3 = df2.dropna(subset=['WT']).sort_values(by=['id','distance'], ascending=True)  # drop NaN value and resorted again.
RV_transect_points_filter_csv = os.path.join(Process_folder,'RV_transect_points_filter.csv')
df3.to_csv(RV_transect_points_filter_csv, index = False)

RV_transect_points_filter = os.path.join(Outputspace, 'RV_transect_points_filtered.shp')  # output file name and location
drv = ogr.GetDriverByName("ESRI Shapefile")  # Set up the shapefile driver
out_ds = drv.CreateDataSource(RV_transect_points_filter)  # Create a data source

prj8 = osr.SpatialReference()  # Create a spatial reference system
prj8.ImportFromWkt(Rasterdata.GetProjectionRef())  # Use same with rasterdata (NAD83(2011) / UTM zone 16N) or ImportFromWkt(Rasterdata.GetProjection())

# Create a layer
out_layer = out_ds.CreateLayer('RV_transect_points_filter', prj8, ogr.wkbPoint)

# reading the points from csv file
with open(RV_transect_points_filter_csv, 'r') as csvfile:
    readerDict = csv.DictReader(csvfile, delimiter=',')
    for field in readerDict.fieldnames:
        new_field = ogr.FieldDefn(field, ogr.OFTString)  # Create a new field with the content of header
        new_field.SetWidth(50)
        out_layer.CreateField(new_field)
    for i, row in enumerate(readerDict):
        # print(row)
        # print(row['x'], row['y'])
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(float(row['x']), float(row['y']))
        feature = ogr.Feature(out_layer.GetLayerDefn())
        feature.SetGeometry(point)
        feature.SetFID(i)
        for field in readerDict.fieldnames:  # Also add x and y coordinates in the shapefile
            j = feature.GetFieldIndex(field)
            feature.SetField(j, row[field])
        out_layer.CreateFeature(feature)

del feature  # destory/delete feature
del out_ds  # Close file

### Step 6 ###########################################################
print("Step 6: Interpolate points to create a raster")
######################################################################
WT_raster_gdal = os.path.join(Outputspace,'WT_raster.tif')
gdal.UseExceptions()
print(extent)

alg_setting = "invdist:power=2.0:smoothing=0.0:radius1=20*Input_Spacing:radius2=20*Input_Spacing:angle=0.0:max_points=12:min_points=2:nodata=0.0"
idw = gdal.Grid(WT_raster_gdal,RV_transect_points_filter, format="GTiff",outputBounds=[extent[0], extent[3], extent[0]+spacing*col_num,extent[3]-spacing*row_num], outputSRS=prj8, width=col_num, height=row_num, outputType=gdal.GDT_Float32, algorithm=alg_setting,zfield='WT')
# caution for the assigned output bounds: [ulx, uly, lrx, lry]
del idw

# Extract IDW data only land regions.
Output_raster = os.path.join(Workspace,'Output_data/WT_raster_extracted.tif')


Input_raster = gdal.Open(WT_raster_gdal,0)
prj = Input_raster.GetProjection()                        # Read projection
# drv = ogr.GetDriverByName("ESRI Shapefile")
dst_ds = ogr.Open(Output_Land_Polygon, 0)  # 0 means read-only. 1 means writeable.
if dst_ds is None:
    sys.exit('Could not open {0}.'.format(dst_ds))
layer = dst_ds.GetLayer()

out_ds = gdal.Warp(Output_raster, WT_raster_gdal,cutlineDSName=Output_Land_Polygon,cropToCutline=True, dstNodata=np.nan)
out_ds.SetProjection(prj)
band = out_ds.GetRasterBand(1)
band.GetMetadata()
out_band = band.SetNoDataValue(NoData_value)  # exclude nodata value
out_band = band.ComputeStatistics(0)    # calculate statistics

del out_ds  # Close file
del Input_raster
del Rasterdata                                                   # Close the raster image

# list of files to be deleted (geopandas products cannot delete now... not sure the reasons.)
# file_list = [file.split(".")[0] + '*' for file in [Baseline_points, Baseline_pre, Offset_line_pre, fn3,RV_transect_points_fishnet, RV_transect_points_fishnet, Smoothed_polyline]]
file_list = [file.split(".")[0] + '*' for file in [Baseline_points, Baseline_pre_points, Offset_line_pre, RV_transect_points_fishnet, RV_transect_points_fishnet, Smoothed_polyline]]

# Delete each file
for pattern in file_list:
    for file in glob.glob(pattern):
        try:
            os.remove(file)
            print(f"{file} has been deleted.")
        except FileNotFoundError:
            print(f"{file} not found.")
        except PermissionError:
            print(f"{file} cannot be deleted due to permission error.")

print("Job Finished ʕ •ᴥ•ʔ")