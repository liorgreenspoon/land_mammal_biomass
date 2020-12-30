import geopandas as gpd
import matplotlib as plt
import numpy as np
import pandas as pd
from matplotlib import pyplot
from matplotlib.colors import LinearSegmentedColormap, hsv_to_rgb
from shapely.geometry import Polygon
import geoplot as gplt
import matplotlib.pyplot
import scipy.interpolate



def get_data(n):
    data = pd.read_csv('new_predictions.csv')
    data = data.drop(['estimated_mass', 'estimated_pop'], axis=1)
    data = data[data.binomial != 'Sus scrofa']       # Wild Boar
    data = data[data.binomial != 'Ursus maritimus']  # Polar bear
    data = data[data.binomial != 'Sus bucculentus']  # EX
    data = data[data.binomial != 'Melomys rubicola'] # EX
    data = data.assign(total_mass=data.AdultBodyMassG * data.pop_density * data.Range,
                       total_mass_density=data.AdultBodyMassG * data.pop_density)
    data = data.sort_values(by='total_mass_density', ascending=False)
    data = data.iloc[0:n - 1]
    geo_data = gpd.read_file('TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp').to_crs("EPSG:6933")
    geo_data = geo_data[geo_data.category != 'EX']
    range_polygons = geo_data.loc[(geo_data['legend'] == 'Extant & Introduced (resident)') |
                                  (geo_data['legend'] == 'Extant & Origin Uncertain (resident)') |
                                  (geo_data['legend'] == 'Extant & Reintroduced (resident)') |
                                  (geo_data['legend'] == 'Extant & Vagrant (seasonality uncertain)') |
                                  (geo_data['legend'] == 'Extant (non breeding)') |
                                  (geo_data['legend'] == 'Extant (resident)') |
                                  (geo_data['legend'] == 'Probably Extant & Origin Uncertain (resident)') |
                                  (geo_data['legend'] == 'Probably Extant (resident)') |
                                  (geo_data['legend'] == 'Reintroduced')]
    range_polygons = range_polygons.merge(data, on='binomial')
    range_polygons = range_polygons.to_crs("EPSG:6933")
    return range_polygons


def get_continent_data():
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world[world.continent != 'Antarctica']
    continents = world.dissolve (by='continent')
    continents = continents.to_crs("EPSG:6933")
    return continents


def get_continent_data_from_file():
    continents = gpd.read_file('map_data/continent_shapefile/continent.shp')
    continents = continents[continents.CONTINENT != 'Antarctica']
    continents = continents.to_crs("EPSG:6933")
    return continents


def get_country_data():
    countries = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    countries = countries[countries.continent != 'Antarctica']
    countries = countries.to_crs("EPSG:6933")
    return countries


def gen_grid(delta):
    outline = get_continent_data()
    xmin, ymin, xmax, ymax = outline.total_bounds
    dx, dy = delta
    xgrid, ygrid = np.meshgrid(np.arange(xmin, xmax, dx), np.arange(ymin, ymax, dy))
    xgrid, ygrid = xgrid.flatten(), ygrid.flatten()
    grid = gpd.GeoDataFrame(geometry=([Polygon([
        [x - dx, y - dy],
        [x - dx, y],
        [x, y],
        [x, y - dy]]) for x, y in zip(xgrid, ygrid)
    ]), crs=outline.crs)

    grid_clip = gpd.clip(grid, outline).reset_index(drop=True)
    grid_clip = grid_clip[grid_clip.geom_type != 'LineString']
    grid_clip = grid_clip[grid_clip.geom_type != 'Point']
    grid_clip = grid_clip[grid_clip.geom_type != 'MultiLineString']
    grid_clip['grid_index'] = grid_clip.index
    grid_clip.grid_index = grid_clip['grid_index'].apply(str)

    return grid_clip


def overlay_and_sum_grid(num_of_species, delta=(500000, 500000)):
    geo_data = get_data(num_of_species)
    # geo_data = geo_data.to_crs("EPSG:6933")
    print("got data")
    grid = gen_grid(delta)
    print("grid generated")
    inter = overlay_and_sum(geo_data, grid)
    print("overlay_and_sum done")
    inter_grouped = inter.groupby(['grid_index']).sum()
    inter_grouped = gpd.GeoDataFrame(inter_grouped.merge(grid, on='grid_index'))
    inter_grouped = inter_grouped.assign(area_km_2=inter_grouped.geometry.area * 10 ** (-6))
    inter_grouped = inter_grouped.assign(total_mass_kg=inter_grouped.total_mass * 10 ** (-3))
    inter_grouped = inter_grouped.assign(total_mass_kg_km2=inter_grouped.total_mass_kg /
                                                           inter_grouped.area_km_2)
    return inter_grouped


def overlay_and_sum_continent(num_of_species):
    geo_data = get_data(num_of_species)
    polygon_layer = get_continent_data_from_file()
    inter = overlay_and_sum(geo_data, polygon_layer)
    inter_grouped = inter.groupby(['CONTINENT']).sum()
    inter_grouped = gpd.GeoDataFrame(inter_grouped.merge(polygon_layer, on='CONTINENT'))
    inter_grouped = inter_grouped.assign(total_mass_Mt=inter_grouped.total_mass * 10 ** (-12))
    return inter_grouped


def overlay_and_sum_country(num_of_species):
    geo_data = get_data(num_of_species)
    polygon_layer = get_country_data()
    inter = overlay_and_sum(geo_data, polygon_layer)
    inter_grouped = inter.groupby(['name']).sum()
    inter_grouped = gpd.GeoDataFrame(inter_grouped.merge(polygon_layer, on='name'))
    inter_grouped = inter_grouped.assign(total_mass_Mt=inter_grouped.total_mass * 10 ** (-12))
    return inter_grouped


def overlay_and_sum(geo_data, layer):
    geo_data = geo_data[~geo_data.is_empty]
    geo_data = geo_data[geo_data.notna()]
    layer = layer[~layer.is_empty]
    layer = layer[layer.notna()]
    inter = gpd.overlay(geo_data, layer, how='intersection')
    inter = inter.assign(inter_area_km_2=inter.geometry.area * 10 ** (-6))
    inter = inter.assign(population=inter.pop_density * inter.inter_area_km_2)
    inter = inter.assign(total_mass=inter.population * inter.AdultBodyMassG)
    return inter


# def gen_custom_cmap():
#     # each row is (z-value, hue, saturation, value/brightness)
#     values = [
#         (0.00, 0.4, 0.0, 1.0),
#         (0.05, 0.4, 0.5, 0.7),
#         (0.10, 0.3, 0.5, 0.6),
#         (0.20, 0.1, 0.5, 0.6),
#         (0.40, 0.0, 0.5, 0.4),
#         (0.60, 0.0, 0.5, 0.25),
#         (0.80, 0.0, 1.0, 0.25),
#         (1.00, 0.0, 1.0, 0.0),
#     ]
#     cdict = {"red": [], "green": [], "blue": []}
#     for z, h, s, v in values:
#         r, g, b = hsv_to_rgb((h, s, v))
#         cdict["red"].append((z, r, r))
#         cdict["green"].append((z, g, g))
#         cdict["blue"].append((z, b, b))
#     cmap = LinearSegmentedColormap('mammals', segmentdata=cdict, N=256)
#     return cmap


def gen_grid_plot(gridded_data, log10=False):
    continents_polygon = gpd.GeoDataFrame(get_continent_data_from_file()).rename(columns={0: "geometry"})
    continents_polygon.crs = "EPSG:6933"
    plt.rcParams["figure.figsize"] = (30, 10)
    font = {'weight': 'normal',
            'size': 40}
    plt.rc('font', **font)
    fig, ax = pyplot.subplots(1, 1)
    ax.axis('off')
    cmap = plt.cm.get_cmap('gist_earth').reversed()
    # cmap = gen_custom_cmap()
    ax.set_label('verbosity coefficient')
    ax.set_title('Wild Mammal Mass Density')
    base = gridded_data.plot(column='total_mass_kg_km2' ,ax=ax, legend=True, cmap=cmap,
                             legend_kwds={'label': r"$kg/km^2$"}, )
    continents_polygon.plot(ax=base, fc='none', ec='grey', linewidth=0.3)


def plot_mass_density_map_from_file():
    gridded_mammal_mass = gpd.read_file('map_data/gridded_mammal_mass/gridded_mammal_mass.shp')
    gridded_mammal_mass = gridded_mammal_mass.rename(columns={"total_ma_1": "total_mass_Mt"})
    gridded_mammal_mass = gridded_mammal_mass.rename(columns={"total_ma_2": "total_mass_kg"})
    gridded_mammal_mass = gridded_mammal_mass.rename(columns={"total_ma_3": "total_mass_kg_km2"})
    gen_grid_plot(gridded_mammal_mass)


def gen_custom_cmap():
    # each row is (z-value, hue, saturation, value/brightness)
    values = [
        (0.00, 0.15, 0.0, 1.0),
        (0.05, 0.15, 0.5, 0.9),
        (0.10, 0.4, 0.5, 0.8),
        (0.15, 0.4, 0.5, 0.7),
        (0.20, 0.2, 0.5, 0.6),
        (0.40, 0.1, 0.5, 0.5),
        (0.80, 0.1, 0.5, 0.4),
        (1.00, 0.0, 1.0, 0.0)
    ]
    cdict = {"red": [], "green": [], "blue": []}
    for z, h, s, v in values:
        r, g, b = hsv_to_rgb((h, s, v))
        cdict["red"].append((z, r, r))
        cdict["green"].append((z, g, g))
        cdict["blue"].append((z, b, b))
    cmap = LinearSegmentedColormap('mammals', segmentdata=cdict, N=256)
    return cmap


