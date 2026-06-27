#!/usr/bin/env python3
"""
NGA_gee_extract.py  —  Earth Engine covariate extraction (LOCAL download).

Pulls environmental covariates for the Nigeria Loa loa high-risk analysis and
writes everything to data/covariates/ on this machine. No Google Drive, no GCS
— only the Earth Engine API (earthengine.googleapis.com) is required.

Outputs:
  data/covariates/covariates_villages.csv   covariates at survey villages (fit)
  data/covariates/covariates_grid.csv        covariates at the 2 km grid (predict)
  data/covariates/maps/<band>.tif            covariate rasters (paper maps)

Covariates (selection happens later in R via horseshoe + projpred):
  elevation, slope (SRTM); LST day, NDVI (MODIS 2016); rainfall total +
  seasonality (CHIRPS 2016); tree cover (Hansen); surface-water occurrence
  (JRC); population (WorldPop 2016).  RAPLOA intentionally excluded.

Auth/project verified working: ee.Initialize(project=EE_PROJECT).
Run with the project venv:  .venv/bin/python code/NGA_gee_extract.py
"""

import os
import io
import zipfile
import urllib.request

import ee
import pandas as pd

EE_PROJECT = "your-earth-engine-project-id"  # set your own GEE Cloud project
YEAR       = 2016
SCALE_M    = 1000
STATES     = ["Abia", "Anambra", "Delta", "Ebonyi", "Imo"]

VILL_CSV   = "data/Emmanuel_cleaned_NGA_data.csv"
GRID_CSV   = "data/prediction_NGA_states_lonlat.csv"   # produced by the R reproject step
OUT_DIR    = "data/covariates"
MAP_DIR    = os.path.join(OUT_DIR, "maps")


def init():
    ee.Initialize(project=EE_PROJECT)


def aoi():
    gaul1 = ee.FeatureCollection("FAO/GAUL/2015/level1")
    return gaul1.filter(ee.Filter.And(
        ee.Filter.eq("ADM0_NAME", "Nigeria"),
        ee.Filter.inList("ADM1_NAME", STATES))).geometry()


def covariate_image():
    d0, d1 = f"{YEAR}-01-01", f"{YEAR}-12-31"

    elev  = ee.Image("USGS/SRTMGL1_003").rename("elevation")
    slope = ee.Terrain.slope(elev).rename("slope")
    lst   = (ee.ImageCollection("MODIS/061/MOD11A2").filterDate(d0, d1)
             .select("LST_Day_1km").mean().multiply(0.02).subtract(273.15)
             .rename("lst_day"))
    ndvi  = (ee.ImageCollection("MODIS/061/MOD13A2").filterDate(d0, d1)
             .select("NDVI").mean().multiply(0.0001).rename("ndvi"))
    rain  = (ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD").filterDate(d0, d1)
             .select("precipitation").sum().rename("rain_total"))

    def monthly(m):
        return (ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD").filterDate(d0, d1)
                .filter(ee.Filter.calendarRange(m, m, "month"))
                .select("precipitation").sum())
    rain_sd = ee.ImageCollection(ee.List.sequence(1, 12).map(monthly)) \
        .reduce(ee.Reducer.stdDev()).rename("rain_seasonality")

    tree  = (ee.Image("UMD/hansen/global_forest_change_2023_v1_11")
             .select("treecover2000").rename("treecover"))
    water = (ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
             .select("occurrence").unmask(0).rename("water_occurrence"))
    pop   = (ee.ImageCollection("WorldPop/GP/100m/pop")
             .filter(ee.Filter.eq("year", YEAR)).mosaic().rename("population"))

    return elev.addBands([slope, lst, ndvi, rain, rain_sd, tree, water, pop]).toFloat()


def sample_points(img, df, out_csv, chunk=1000):
    """Sample img at df rows (long/lat) in chunks, write a tidy CSV locally."""
    bands = img.bandNames().getInfo()
    rows = []
    for start in range(0, len(df), chunk):
        sub = df.iloc[start:start + chunk]
        feats = [ee.Feature(ee.Geometry.Point([float(r.long), float(r.lat)]),
                            {"row_id": int(i)}) for i, r in sub.iterrows()]
        fc = img.reduceRegions(ee.FeatureCollection(feats),
                               ee.Reducer.first(), SCALE_M)
        for f in fc.getInfo()["features"]:
            p = f["properties"]
            rows.append({"row_id": p["row_id"], **{b: p.get(b) for b in bands}})
        print(f"  sampled {min(start + chunk, len(df))}/{len(df)}")
    cov = pd.DataFrame(rows).set_index("row_id")
    out = df.join(cov, how="left")
    out.to_csv(out_csv, index=False)
    print("wrote", out_csv, out.shape)


def download_maps(img, region):
    os.makedirs(MAP_DIR, exist_ok=True)
    for band in img.bandNames().getInfo():
        url = img.select(band).clip(region).getDownloadURL(
            {"region": region, "scale": SCALE_M, "format": "GEO_TIFF",
             "crs": "EPSG:4326"})
        dest = os.path.join(MAP_DIR, f"{band}.tif")
        try:
            urllib.request.urlretrieve(url, dest)
            print("map:", dest)
        except Exception as e:                       # some servers return a zip
            data = urllib.request.urlopen(url).read()
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                z.extractall(MAP_DIR)
            print("map (zip):", band, "->", MAP_DIR)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    init()
    region = aoi()
    img = covariate_image()

    vill = pd.read_csv(VILL_CSV)
    vill = (vill[vill["NAME_1"].isin(STATES)]
            .dropna(subset=["long", "lat"])
            .drop_duplicates(["NAME_1", "NAME_2", "Vill"])
            .reset_index(drop=True))
    grid = pd.read_csv(GRID_CSV)

    print(f"villages: {len(vill)}  grid: {len(grid)}")
    sample_points(img, vill, os.path.join(OUT_DIR, "covariates_villages.csv"))
    sample_points(img, grid, os.path.join(OUT_DIR, "covariates_grid.csv"))
    download_maps(img, region)
    print("Done. Covariates in", OUT_DIR)


if __name__ == "__main__":
    main()
