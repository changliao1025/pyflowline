###########################
Frequently Asked Questions
###########################

1. Why my `conda` cannot create environment?
   
   Turn off the VPN or bypass it.

2. Why import `GDAL` failed?
   
   Consider using the `conda-forge` channel.

3. `proj` related issue https://github.com/OSGeo/gdal/issues/1546, 
   
   Make sure you correctly set up the `PROJ_LIB`

   Because the `GDAL` library is used by this project and the `proj` library is often not configured correctly automatically. 
   On Linux or Mac, you can set it up using the `.bash_profile` such as:

   Anaconda:

   `export PROJ_LIB=/people/user/.conda/envs/hexwatershed/share/proj`
   `export PROJ_LIB=$HOME/opt/anaconda3/envs/pyflowline/share/proj`

   Miniconda:

   `export PROJ_LIB=/opt/miniconda3/envs/hexwatershed/share/proj`

4. I am getting errors using the plot functions

   If you receive an error related to `GeoAxesSubplot`, make sure you have cartopy version 0.21.0 installed in your environment. Optionally, downgrade matplotlib to 3.5.2. 

5. What if my model doesn't produce the correct or expected answer?
   
   Answer: There are several hidden assumptions within the workflow. For example, if you provide the DEM and river network for two different regions, the program won't be able to tell you that. A visual inspection of your data is important.
   
   Optionally, you can turn on the `iFlag_debug` option in the configuration file to output the `intermediate files`.