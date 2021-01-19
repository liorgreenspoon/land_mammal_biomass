import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon
from
תוimport zonal_stats

path = 'shp_files/'
annual_mean_temp_str = 'data/annual_mean_temp_projected.tif'
prec_warmest_quarter_str = 'data/prec_warmest_quarter_projected.tif'
temp_seasonality_str = 'data/temp_seasonality_projected.tif'
npp_str = 'data/npp_projected.tif'
prec_seasonality_str = 'data/prec_seasonality_projected.tif'
stat = "mean"
spatial_res_grid = 20000
spatial_res = 100*10**6


def polygon_to_grid(delta, polygon):
    xmin, ymin, xmax, ymax = polygon.bounds
    dx, dy = delta
    xgrid, ygrid = np.meshgrid(np.arange(xmin, xmax, dx), np.arange(ymin, ymax, dy))
    xgrid, ygrid = xgrid.flatten(), ygrid.flatten()
    grid = gpd.GeoDataFrame(geometry=([Polygon([
        [x - dx, y - dy],
        [x - dx, y],
        [x, y],
        [x, y - dy]]) for x, y in zip(xgrid, ygrid)
    ]), crs='epsg:6933')
    grid_clip = gpd.clip(grid, polygon).reset_index(drop=True)
    grid_clip = grid_clip[grid_clip.geom_type != 'LineString']
    grid_clip = grid_clip[grid_clip.geom_type != 'Point']
    grid_clip = grid_clip[grid_clip.geom_type != 'MultiLineString']
    grid_clip['grid_index'] = grid_clip.index
    grid_clip.grid_index = grid_clip['grid_index'].apply(str)
    grid_clip = grid_clip.assign(area=grid_clip.area * 10 ** (-6))

    return grid_clip



def get_gridded_data(one_species):
    binomial = one_species.name.replace(' ', '_')
    polygon = one_species.geometry
    if polygon.area > spatial_res:
        grid = polygon_to_grid((spatial_res_grid, spatial_res_grid), polygon)
    else:
        grid = one_species
    grid.to_file(path + binomial + '.shp')
    grid = grid.assign(binomial=binomial,
                       annual_mean_temp=pd.DataFrame(zonal_stats(path + binomial + '.shp',
                                                                 annual_mean_temp_str,
                                                                 stats=stat, all_touched=True)),
                       prec_warmest_quarter=pd.DataFrame(zonal_stats(path + binomial + '.shp',
                                                                     prec_warmest_quarter_str,
                                                                     stats=stat, all_touched=True)),
                       temp_seasonality=pd.DataFrame(zonal_stats(path + binomial + '.shp',
                                                                 temp_seasonality_str,
                                                                 stats=stat, all_touched=True)),
                       npp=pd.DataFrame(zonal_stats(path + binomial + '.shp',
                                                    npp_str,
                                                    stats=stat, all_touched=True)),
                       prec_seasonality=pd.DataFrame(zonal_stats(path + binomial + '.shp',
                                                                 prec_seasonality_str,
                                                                 stats=stat, all_touched=True))
                       )
    grid = grid[['grid_index', 'area', 'annual_mean_temp',
                 'prec_warmest_quarter', 'temp_seasonality', 'prec_seasonality', 'npp', 'binomial']].dropna()
    grid.to_csv('results/' + binomial + '.csv')
    print(binomial + ' done')
