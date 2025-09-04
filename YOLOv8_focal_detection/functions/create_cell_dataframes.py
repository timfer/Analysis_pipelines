import os
import pickle
import pandas as pd
import re
import numpy as np

from utils.dir_paths import paths
from utils.variables import dropped_planes
from functions.dir_utils import pattern

def create_template_df():
    '''
    Input:
    - None required.

    Action:
    - Creates an empty pandas DataFrame with predefined columns to store information about focal cells. 
    - These columns include identifiers like 'Full_Cell_ID', experimental conditions, and cell count, among others.
    - The structure of this DataFrame is designed to facilitate easy addition and manipulation of cell data across different experimental conditions and imaging planes.

    Output:
    - An empty DataFrame with columns for storing cell data.
    '''
    df = pd.DataFrame({
        'Full_Cell_ID': pd.Series(dtype='str'),
        'exp_ID': pd.Series(dtype='str'),
        'plate_code': pd.Series(dtype='str'),
        'well_code': pd.Series(dtype='str'),
        'cell_code': pd.Series(dtype='str'),
        'Cluster_names': pd.Series(dtype='str'),
        'Cell_subclass': pd.Series(dtype='str'),
        'Condition': pd.Series(dtype='str'),
        'nCount_RNA': pd.Series(dtype='int'),
        'nFeature_RNA': pd.Series(dtype='str'),
        'centroids': pd.Series(dtype='object'),
        'channels': pd.Series(dtype='object'),
        'cell_count': pd.Series(dtype='float64')
    })
    return df



def split_full_cell_id(df, column_name='Full_Cell_ID'):
    '''
    Input:
    - df: A pandas DataFrame containing a column 'Full_Cell_ID' that has compound identifiers for cells.

    Action:
    - Splits the 'Full_Cell_ID' column into multiple new columns ('exp_ID', 'plate_code', 'well_code', 'cell_code') for detailed identification.
    - Each part of the 'Full_Cell_ID' is expected to be separated by underscores and is parsed accordingly.
    - New columns for 'centroids' and 'channels' are initialized as NaN, preparing the DataFrame for further data assignment.

    Output:
    - The modified DataFrame with additional columns derived from 'Full_Cell_ID' and placeholders for cell data.
    '''
    df['exp_ID'] = df[column_name].str.extract(r'(?<![A-Z])(?!CC)([A-Z]{2}\d{2,3})(?![A-Z])') #Will extract the pattern of the experimental code
    df['plate_code'] = df[column_name].str.extract(r'((?<=_)P\d{1,3})') #Will extract the platecode
    df['cell_code'] = df[column_name].str.extract(r'((?<=_)CC\d{1,3})') #Extracts the cellcode
    df['well_code'] = df[column_name].str.extract(r'((?<=_)WC\d{1,3})') #Extracts the wellcode
    for new_col in ['centroids', 'channels']:
        df[new_col] = np.nan
    return df

def create_cell_dataframe():
    '''
    Input:
    - Assumes the existence of a global variable `centroid_dico` containing centroid data, and paths `resultPath` and `yoloPath`.

    Action:
    - Loads a dictionary of centroid data from a pickle file located at `resultPath`.
    - Creates a template DataFrame to structure the cell data for analysis.
    - Iterates through the `yoloPath` directory, constructing identifiers for cells based on directory names and files, and populating the DataFrame with these identifiers.
    - Processes the loaded centroid data, splitting identifiers into detailed components, and filling in cell-specific data (centroids, channels, cell counts).
    - Differentiates cells based on their count into single, none, and multi-cell categories, saving separate CSV files for each category for easy access and analysis.

    Output:
    - A detailed pandas DataFrame `df` containing all processed cell data, alongside saving specific subsets of the data as CSV files categorized by cell count.
    '''
    #Load in the centroid dictionary established and saved in the resultPath
    with open(os.path.join(paths['resultPath'], 'centroid_dico.pkl'), 'rb') as pkl_file:
        centroid_dico = pickle.load(pkl_file)

    #Transform the centroid dictionary into a dataframe
    df = create_template_df()

    exp_ID = re.search(r'(?<![A-Z])(?!CC)([A-Z]{2}\d{2,3})(?![A-Z])', (paths['yoloPath'].replace('\\', '/')).split('/')[-2]).group(0)
    new_rows = []

    #Now since this is an empty dataframe, we need to iterate through our different image files to create a 'Full_Cell_ID'
    for pc_cc in os.listdir(paths['yoloPath']):
        if not re.search(pattern, pc_cc): #Same pattern as that established in dir_utils
            continue

        plate_code, cell_code = str(pc_cc.split('_')[0]), str(pc_cc.split('_')[1])
        if os.path.isdir(os.path.join(paths['yoloPath'], pc_cc)):
            for cell_id in os.listdir(os.path.join(paths['yoloPath'], pc_cc)):
                wc = f"WC{str(int(cell_id.split('_')[-1])+1)}"
                full_cell_id = f'{exp_ID}_{plate_code}_{wc}_{cell_code}'
                # Append a new row to the dataframe
                new_rows.append({'Full_Cell_ID': full_cell_id})

    if new_rows:
        df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    #First we split the 'Full_Cell_ID' strings into 'experiment_ID', 'plate_code', 'well_code', 'cell_code'
    df = split_full_cell_id(df)

    for pc_cc, cell_data in centroid_dico.items():
        for cell_id, centroid in cell_data.items():
            plate_code, cell_code = str(pc_cc.split('_')[0]), str(pc_cc.split('_')[1])
            cell_count = centroid['cell_count']
            wc = f"WC{str(int(cell_id.split('_')[-1])+1)}" #We have to add 1 because in this case WC are 1 indexed
            condition = (df['plate_code'] == plate_code) & (df['cell_code'] == cell_code) & (df['well_code'] == wc)
         
            if cell_count == None or not centroid['data']:
                print('No cells')
                # If cell_count is 0 or centroid['data'] is empty, skip this iteration
                df.loc[condition, 'cell_count'] = 0
                continue
            centroids = []
            channels = []
            for _, data in centroid['data'].items():
                print(data)
                #Sometimes "fake cells" get by and we need to correct for this, here we check if centroid['data'] is empty
                #If the focal channel is 'ch1', this is highly likely to be trash/background that is mislabeled as 
                if data['centroid'][0].any():
                    if data['channel'][0] not in dropped_planes:
                        centroids.append(tuple(data['centroid'][0]))
                        #print(f'idx: {idx}, {data[1]}')
                        channels.append(str(data['channel'][0]))
            
            df.loc[df.index[condition], 'centroids'] = df.loc[df.index[condition], 'centroids'].apply(lambda x: centroids if centroids else None)
            df.loc[df.index[condition], 'channels'] = df.loc[df.index[condition], 'channels'].apply(lambda x: channels if channels else None)
            df.loc[condition, 'cell_count'] = len(centroids)
    #Save dataframe for further use
    df.to_csv(os.path.join(paths['resultPath'], 'complete_df.csv'), index=False)

    single_cell_df = df[df['cell_count'] == 1]
    single_cell_df.to_csv(os.path.join(paths['resultPath'], 'single_cell_df.csv'), index=False)

    no_cell_df = df[df['cell_count'] == 0]
    no_cell_df.to_csv(os.path.join(paths['resultPath'], 'no_cell_df.csv'), index=False)

    multi_cell_df = df[df['cell_count'] > 1]
    multi_cell_df.to_csv(os.path.join(paths['resultPath'], 'multi_cell_df.csv'), index=False)

    return df