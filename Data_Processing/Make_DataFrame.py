import glob
import pandas as pd
import os
from Data_Processing.EAG_DataProcessing_Library import open_wave
#dir='Manduca2Reciever_Data/Buttered/Sum/Pipette/'
#iles = glob.glob(f"{dir}/*/*.csv")
#or file in files:
 ##   meta_data = os.path.basename(file.lower()).split("_")
  #  print(meta_data[7][:5])

def EAG_Cov_df_build(DIR):
    """
    Builds a Pandas DataFrame from all CSV files in a directory.

    Args:
        DIR (str): Path to the directory containing the CSV files.

    Returns:
        A Pandas DataFrame containing the wave data from all CSV files in the directory.
        The DataFrame has one row per file, with columns for the wave data, label, concentration, and date.
    """
    files = glob.glob(f"{DIR}/*/*.csv")
    df_list = []  # to store each DataFrame and concatenate at the end
    for file in files:
        if '_p_' in file.lower():
            # Open the file
            wave_data = open_wave(file)

            # Convert wave_data to a DataFrame
            df_wave = pd.DataFrame(wave_data).T

            # Get metadata from the file name
            meta_data = os.path.basename(file.lower()).replace('.csv', '').split("_")
            if meta_data[2]== 'p':
                meta_data[2]='pipette'

            # Prepare a DataFrame for meta_data
            df_meta = pd.DataFrame({
                'Wave': [meta_data[7][:5]],
                'Record_Channel': [meta_data[6]],
                'Delivery': [meta_data[2]],
                'Line': [meta_data[1]],
                'Label': [meta_data[4]],
                'Concentration': [meta_data[3]],
                'Date': [meta_data[0]],
            })

            # Concatenate df_wave and df_meta along columns
            df_file = pd.concat([df_wave, df_meta], axis=1)

            # Change the index to be the filename without '.csv'
            df_file.index = [os.path.basename(file).replace('.csv', '')]

            # Append the DataFrame to our list
            df_list.append(df_file)

    # Concatenate all DataFrames in df_list to create the final DataFrame
    df = pd.concat(df_list)

    return df




dir= '../Manduca2Reciever_Data/Buttered/Norm_Independent/Pipette/'
df = EAG_Cov_df_build(dir)
#N_Full.loc[N_Full['label']  == str('healthy1k'), 'label'] = 'Healthy'
df.loc[df['Label'] == str('linelool'), 'Label'] = 'linalool'
df.to_csv('Manduca2Reciever_Data/DataFrames/NoQC/070723Buttered_Norm_Independent_Pipette.csv')