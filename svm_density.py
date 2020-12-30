import pandas as pd
# from tqdm.auto import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
import numpy as np
import math
from sklearn.metrics import mean_squared_error

# tqdm.pandas()
pd.set_option('display.max_columns', 10000)


def read_data() -> pd.DataFrame:
    data_file_names = dict(
        mammal_labeled='svr_data/new_svr_data.csv',
        mammal_not_labeled='svr_data/Mammal_AI_not_global_no_overlap.csv',
        clades_labeled='svr_data/clade_df_labeled.csv',
        clades_not_labeled='svr_data/clade_df.csv'
    )
    df = pd.concat(
        [pd.read_csv(data_file_names['mammal_labeled']),
         pd.read_csv(data_file_names['mammal_not_labeled'])],
        keys=['labeled', 'not_labeled'], names=['labeled'], sort=True
    ).set_index('Binomial', append=True).rename_axis(['labeled', 'id', 'binomial'])
    common_names = pd.read_csv('svr_data/common_names.csv')
    common_names = common_names.set_index('binomial').drop(['Unnamed: 0'], axis=1)
    df = df.join(common_names, on='binomial')
    clades = pd.concat(
        [pd.read_csv(data_file_names['clades_labeled']),
         pd.read_csv(data_file_names['clades_not_labeled'])]
    ).drop(['clade1', 'clade2', 'clade4', 'clade7', 'clade8', 'clade9'], axis=1)\
        .drop_duplicates().set_index('binomial')
    df = df.join(clades, on='binomial')

    return df


def preproc_data(data_raw: pd.DataFrame):
    numerical_columns_rename = dict(Range='range', AdultBodyMassG='body_mass', GenerationLengthD='gen_length',
                                    total_pop='population', pop_density='density')
    data = (data_raw.rename(columns=numerical_columns_rename))
    data['RedListStatus'] = data['RedListStatus'].replace(['DD', 'NotAssigned'], 'LC')
    data = data[data.RedListStatus != 'EX']
    data = data.assign(log_range=np.log10(data.range),
                       log_population=np.log10(data.population),
                       log_body_mass=np.log10(data.body_mass),
                       log_density=np.log10(data.density))
    red_list_status = pd.get_dummies(data.RedListStatus, prefix='red_list_status')
    clades = pd.get_dummies(data.clade, prefix='clade')
    data = pd.concat([data, red_list_status, clades], axis=1)
    data = data.drop(['Family', 'Genus', 'Order', 'Range_m', 'RedListStatus', 'TrophicLevel',
                      'population', 'gen_length', 'label', 'population'], axis=1)
    cat_features = list(red_list_status.columns)+list(clades.columns)
    cont_features = ['log_range', 'log_body_mass']
    label_name = ['log_density']
    return [data, cat_features,  cont_features,label_name]


def scale_cont(dataset, cont_features, label_name):
    feature_scaler = StandardScaler()
    label_scaler = StandardScaler()
    scaled_features = feature_scaler.fit_transform(dataset[cont_features].values)
    scaled_features = pd.DataFrame(scaled_features, index=dataset.index, columns=cont_features)
    scaled_labels = label_scaler.fit_transform(dataset[label_name].values)
    scaled_labels = pd.DataFrame(scaled_labels, index=dataset.index, columns=label_name)
    return [feature_scaler, label_scaler, scaled_features, scaled_labels]


def svr(labeled_data, cat_features, cont_features, label_name, random_state):
    feature_scaler, label_scaler, scaled_features, scaled_labels = \
        scale_cont(labeled_data, cont_features, label_name)
    scaled_df = pd.concat([scaled_features, labeled_data[cat_features]], axis=1)
    x_train, x_test, y_train, y_test = train_test_split(scaled_df, scaled_labels, test_size=0.1,
                                                        random_state=random_state)
    regressor = SVR(kernel='rbf', epsilon=0.1)
    regressor.fit(x_train, y_train.values.ravel())
    predictions_df = pd.DataFrame(label_scaler.inverse_transform(regressor.predict(x_test)),
                                  index=x_test.index, columns=['predictions'])
    test_set = pd.DataFrame(x_test.apply(lambda row: row.name, axis=1))
    result_df = test_set.join(labeled_data).drop([0], axis=1)
    result_df = result_df.join(predictions_df)
    if label_name[0] == 'log_density':
        rmse = math.sqrt(mean_squared_error(result_df.predictions, result_df[label_name]))
    elif label_name[0] == 'log_population':
        rmse = math.sqrt(mean_squared_error(result_df.predictions-result_df.log_range,
                                            result_df.log_population-result_df.log_range))
    return [result_df, rmse]

def update_result_df(result_df):
    result_df['prediction_mean'] = result_df.loc[:, result_df.columns == 'predictions'].mean(numeric_only=True, axis=1)
    result_df = result_df.drop(['predictions'], axis=1)
    return result_df


def svr_predict(data, cat_features, cont_features, label_name):
    feature_scaler, label_scaler, scaled_features, scaled_labels = \
        scale_cont(data, cont_features, label_name)
    scaled_df = pd.concat([scaled_features, data[cat_features]], axis=1)
    regressor = SVR(kernel='rbf', epsilon=0.1)
    regressor.fit(scaled_df.loc['labeled'], scaled_labels.loc['labeled'].values.ravel())
    log_predictions = label_scaler.inverse_transform(regressor.predict(scaled_df.loc['not_labeled']))
    predictions = np.power(10, (log_predictions + (0.69**2) / 2 * np.log(10)))     # E(x)= 10[μ+22*ln(10)]
    predictions_df = pd.DataFrame(predictions, index=scaled_df.loc['not_labeled'].index, columns=['predictions'])
    predictions_df = data.loc['not_labeled'].join(predictions_df)
    if label_name[0] =='log_density':
        predictions_df = predictions_df.assign(population=predictions_df.predictions *
                                                             predictions_df.range)
        predictions_df = predictions_df.assign(total_mass_Mt=predictions_df.population *
                                                             predictions_df.body_mass * 10 ** (-12))
        # print('DENSITY')
    elif label_name[0] =='log_population':
        predictions_df = predictions_df.assign(total_mass_Mt=predictions_df.predictions *
                                                             predictions_df.body_mass * 10**(-12))
        # print('POP')
    return predictions_df

