import pandas as pd


def subset_and_rename(df, col_name, common_names_dict, df_type):
    if dict != 0:
        df = df.replace({col_name: common_names_dict})
    df = df.rename({col_name: 'Level_2', 'total_mass': 'mass'}, axis=1)
    df = df.groupby('Level_2').agg('sum').reset_index()
    df['Level_1'] = df_type
    df = df[['Level_1', 'Level_2', 'mass']]
    return df


def main():
    terrestrial = pd.read_csv('voronoi_diagram_data/Terrestrial.csv')
    marine = pd.read_csv('voronoi_diagram_data/Marine.csv')
    domesticated = pd.read_csv('voronoi_diagram_data/Domesticated_and_human.csv')

    marine_dict = {'Balaenidae': 'Baleen whales',
                   'Balaenopteridae': 'Baleen whales',
                   'Neobalaenidae': 'Baleen whales',
                   'Eschrichtiidae': 'Baleen whales',
                   'Physeteridae': 'Sperm whales',
                   'Kogiidae': 'Sperm whales',
                   'Ziphiidae': 'Beaked whales',
                   'Monodontidae': 'Dolphins and Porpoises',
                   'Platanistidae': 'Dolphins and Porpoises',
                   'Iniidae': 'Dolphins and Porpoises',
                   'Lipotidae': 'Dolphins and Porpoises',
                   'Delphinidae': 'Dolphins and Porpoises',
                   'Pontoporiidae': 'Dolphins and Porpoises',
                   'Phocoenidae': 'Dolphins and Porpoises',
                   'Trichechidae': 'Sea cows',
                   'Dugongidae': 'Sea cows',
                   'Otariidae': 'Otaries',
                   'Phocidae': 'True seals',
                   'Ursidae': 'Other Carnivores',
                   'Mustelidae': 'Other Carnivores',
                   'Odobenidae': 'Other Carnivores'
                   }

    terrestrial_dict = {'Proboscidae': 'Elephants',
                        'Perissodactyla': 'Odd-toed hoofed mammals',
                        'Diprotodontia': 'Pouched mammals',
                        'Cetartiodactyla': 'Even-toed hoofed mammals',
                        'Carnivora': 'Carnivores',
                        'Primates': 'Primates',
                        'Eulipotyphla': 'Insect-eaters',
                        'Chiroptera': 'Bats',
                        'Rodentia': 'Rodents'
                        }

    marine_dict = {k.upper(): v for k, v in marine_dict.items()}
    terrestrial_dict = {k.upper(): v for k, v in terrestrial_dict.items()}

    # All terrestrial Orders not in dictionary are defined as "Other"
    terrestrial.loc[~terrestrial['Order'].isin(terrestrial_dict.keys()), 'Order'] = 'Other'

    marine = subset_and_rename(marine, 'Family', marine_dict, 'Marine')
    terrestrial = subset_and_rename(terrestrial, 'Order', terrestrial_dict, 'Terrestrial')
    domesticated = subset_and_rename(domesticated, 'name', 0, 'Domesticated')

    all_species = pd.concat([marine, terrestrial, domesticated])
    all_species.loc[all_species['Level_2'] == 'Human', 'Level_1'] = 'Human'
    all_species.to_csv('hierarchy_for_voronoi.csv')


if __name__ == "__main__":
    main()




