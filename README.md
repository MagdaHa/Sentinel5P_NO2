# Sentinel-5P: Processing and handling NO2 measured values in R and ArcGIS

Sentinel-5P is a satellite that is operating within ESA's Copernicus program since 2017. The goal of the mission is a very dense scheduled operational monitoring of the atmosphere. Using the TROPOMI instrument on board, various parameters are measured, such as tropospheric pollutants, greenhouse gases, aerosols, and cloud cover. NO2 concentrations are also provided as a separate product, which can be downloaded as NetCDF.

![S5P product](https://github.com/MagdaHa/Sentinel5P_NO2/blob/master/img/S5P_product.png)

More information about this product can be found in the [ESA user manual](https://sentinel.esa.int/documents/247904/2474726/Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide).

## R script
With the help of the R package [GetSpatialData](https://github.com/16EAGLE/getSpatialData) various satellite data available on the Copernicus Open Access Hub can be downloaded automatically. Conditions like recording period, satellite, product type and area of interest can be defined.

NO2 data is stored in the Sentinel-5P product type L2__NO2___, therefore only this data is of interest and should be downloaded for the desired period.
The subsequently downloaded NetCDF files contain several variables, such as longitude and latitude, different viewing angels and various nitrogen dioxide related variables. Product variable 6 (nitrogendioxide_tropospheric_column) displays the NO2 concentration in the troposphere, which is used for the following analysis.

Using the [S5Processor](https://github.com/MBalthasar/S5Processor) package, this variable can be queried and saved as GeoTiff. In addition, cloud masking is applied to prevent outliers. Product variable qa_value indicates the error probability of the measured values (e.g. due to cloud cover or snow). By blending the two variables and the resulting elimination of errors (qa_value < 0.5 or 0.75 = NoData) a corrected nitrogen dioxide layer is obtained (unit molecules/cm²).

![R workflow](https://github.com/MagdaHa/Sentinel5P_NO2/blob/master/img/R_workflow.png)

#### How are the result layers supposed to be interpreted?
This script focuses on the technology used and not on the scientific interpretation of the results. Nevertheless:
* Values provide information about NO2 concentrations in the troposphere (up to 15 km above the earth's surface), but not at a specific location near the ground (e.g. a busy         road)
* NO2 concentrations are based on different influencing factors: Winds influence the airflow of NO2 from or in neighboring countries; cloud cover falsifies measured values and       must be masked; low temperatures prevent the air masses from rising and thus also the rise of NO2; thunderstorms: lightning as a source of nitrogen oxides
* Averaging over a certain period bypasses NoData values and eliminates meteorological effects that occur selectively
  
  
A detailed StoryMap of Esri Germany with illustrations and the Sentinel-5P image service can be found [here](https://storymaps.arcgis.com/stories/ffb2678bf09f466b9744d30c5fb902a2) (only in German).


## ArcGIS Pro Toolbox
The toolbox integrate the steps of this workflow using R Bridge and native arcpy. The toolbox does not do the data aggregation but relies on the ArcGIS Pro capabilities of flexible aggregation (for days, weeks, months, …) through the Multidimensional Toolset in ArcGIS Pro.


To make sure both the R script and the toolbox run without problems, the following packages must be installed in R: devtools, getSpatialData, raster, sf, sp, dismo, geosphere, rgeos, S5Processor, ncdf4, ggplot2, maptools, rgdal, lubridate, cat. 
For the toolbox, a working version of the [R-bridge](https://github.com/R-ArcGIS/r-bridge-install) needs to be installed, too.

