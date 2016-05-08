# NCEP REANALYSIS 2

This folder contains _NCEP-DOE AMIP-II Reanalysis (AKA Reanalysis 2)_ files to
research. Files are downloaded from [NCEP Reanalysis 2](http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html).

To download files below:
run sahra.Rproj in rstudio, source filehelper.r and issue the command

```R
download_ncep_R2()
```

###### Common Properties:
* **Latitude :** 0 - 60N
* **Longitude :** 10E - 55E
* **Spatial Resolution :** 2.5x2.5 degree for pressure levels.
* **Pressure Levels :** 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200
* **Time :** 2008-01-01 00:00:00 - 2012-12-31 18:00:00

###### Files:

* **r2-pres-4-hgt.nc :** Geopotential height at pressure levels / 4 times daily.
* **r2-pres-4-omega.nc :** Omega (dP/dt) at pressure levels / 4 times daily.
* **r2-pres-4-u.nc :** u component of the wind at pressure levels / 4 times daily.
* **r2-pres-4-v.nc :** v component of the wind at pressure levels / 4 times daily.
* **r2-pres-daily-hgt.nc :** Geopotential height at pressure levels / Daily mean.
* **r2-pres-daily-omega.nc :** Omega (dP/dt) at pressure levels / Daily mean.
* **r2-pres-daily-u.nc :** u component of the wind at pressure levels / Daily mean.
* **r2-pres-daily-v.nc :** u component of the wind at pressure levels / Daily mean.
* **r2-surf-4-mslp.nc :** MSLP at surface / 4 times daily.
* **r2-surf-4-u.nc :** u component of the wind at 10m / 4 times daily.
* **r2-surf-4-v.nc :** v component of the wind at 10m / 4 times daily.
* **r2-surf-daily-mslp.nc :** MSLP at surface / Daily mean.
* **r2-surf-daily-u.nc :** u component of the wind at 10m / Daily mean.
* **r2-surf-daily-v.nc :** u component of the wind at 10m / Daily mean.
