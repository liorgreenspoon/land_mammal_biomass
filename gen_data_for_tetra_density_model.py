import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon
from rasterstats import zonal_stats


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
    path = 'tetra_density_data/shp_files/'
    stat = "mean"
    spatial_res_grid = 20000
    spatial_res = 100*10**6

    annual_mean_temp_str = 'tetra_density_data/annual_mean_temp_projected.tif'
    prec_warmest_quarter_str = 'tetra_density_data/prec_warmest_quarter_projected.tif'
    temp_seasonality_str = 'tetra_density_data/temp_seasonality_projected.tif'
    npp_str = 'tetra_density_data/npp_projected.tif'
    prec_seasonality_str = 'tetra_density_data/prec_seasonality_projected.tif'
    binomial = one_species.binomial.replace(' ', '_')
    polygon = one_species.geometry
    if polygon.area > spatial_res:
        grid = polygon_to_grid((spatial_res_grid, spatial_res_grid), polygon)
    else:
        grid = one_species
    try:
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
        grid = grid[['grid_index', 'area', 'annual_mean_temp', 'prec_warmest_quarter', 'temp_seasonality',
                     'prec_seasonality', 'npp', 'binomial']].dropna()
        pd.DataFrame(grid).to_csv('tetra_density_data/results/' + binomial + '.csv')
        print(binomial + ' done')
    except:
        print("Couldn't process "+ binomial)


def main():
    mammals_df_geo = gpd.read_file('TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp').to_crs("EPSG:6933")
    ## Uncomment the section below to run for just a subset
    # subset = pd.read_csv('missing_tetra_density_species.csv')
    # mammals_df_geo = gpd.GeoDataFrame(subset.merge(mammals_df_geo, on='binomial', how='left'))
    mammals_df_geo = mammals_df_geo.reset_index()
    mammals_df_geo = mammals_df_geo.dissolve(by='binomial').reset_index()
    mammals_df_geo.apply(get_gridded_data, axis=1)


if __name__ == "__main__":
    main()
