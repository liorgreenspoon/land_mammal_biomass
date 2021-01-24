import pandas as pd
import numpy as np
global santini_models
santini_models = pd.read_csv('tetra_density_data/santini_models.csv').fillna(0)


def calculate_density(x):
    santini_models_copy = santini_models.copy()
    santini_models_copy = santini_models_copy.multiply([1, x.log_BM, x.log_BM_2, x.log_BM_3, x.Herbivore,
                                                        x.Omnivore, x.log_npp, x.log_npp_2, x.Pwq, x.Pwq_2,
                                                        x.Pcv, x.Pcv_2],
                                                       axis='columns')
    santini_models_copy['density'] = 10**santini_models_copy.sum(axis=1)
    return np.mean(santini_models_copy['density'])


def main():
    env_params = pd.read_csv('tetra_density_data/env_params_for_tetra_density.csv')
    env_params.binomial = env_params.binomial.str.replace('_', ' ')
    env_params.loc[env_params.npp < 0, 'npp'] = 1
    env_params['log_npp'] = np.log10(env_params['npp'])
    env_params = env_params[env_params.annual_mean_temp >= -100]

    species_params = pd.read_csv('tetra_density_data/all_species_list.csv')
    species_params = pd.concat([species_params, (pd.get_dummies(species_params.TrophicLevel))], axis=1)
    species_params = species_params.drop(['Range', 'Red List status', 'Order', 'TrophicLevel',
                                          'NotAssigned', 'Carnivore', 'Insectivore'], axis=1)
    species_params['log_BM'] = np.log10(species_params['AdultBodyMass_g'])

    all_params = env_params.merge(species_params, on='binomial').drop(['Unnamed: 0'], axis=1)
    all_params = all_params.rename({'annual_mean_temp': 'Tm', 'prec_seasonality': 'Pcv',
                                    'prec_warmest_quarter': 'Pwq', 'area': 'area_km_2'}, axis=1)
    all_params['area_ha'] = all_params['area_km_2'] / 100
    all_params = all_params.assign(log_BM_2=all_params.log_BM ** 2, log_BM_3=all_params.log_BM ** 3,
                                   Pwq_2=all_params.Pwq ** 2,
                                   log_npp_2=all_params.log_npp ** 2,
                                   Pcv_2=all_params.Pcv ** 2)

    all_params['density'] = all_params[['log_BM', 'log_BM_2', 'log_BM_3', 'Herbivore', 'Omnivore', 'log_npp',
                                        'log_npp_2', 'Pwq', 'Pwq_2', 'Pcv', 'Pcv_2']].apply(calculate_density, axis=1)
    all_params['santini_population'] = all_params['density'] * all_params['area_ha'] * 10 ** (1.0 ** 2 / 2 * np.log(10))
    all_params['santini_total_mass_Mt'] = all_params['santini_population'] * all_params['AdultBodyMass_g'] * 10 ** (-12)

    output_df = all_params[['binomial', 'santini_population', 'santini_total_mass_g']]
    output_df.groupby(['binomial']).sum().to_csv('santini_output.csv')


if __name__ == "__main__":
    main()
