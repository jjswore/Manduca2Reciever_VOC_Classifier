from Data_Processing.EAG_DataProcessing_Library import *
import matplotlib.pyplot as plt
import csv

def get_EAGs(data):
    """
        Extracts the timeseries data for each channel in a given dataset.

        Args:
            data (list): A 2D list representing the dataset containing EAG data.
                         The data should be organized such that each row represents
                         a time point, and each column represents a different channel.
                         The first row contains the channel names, and the following
                         rows contain the corresponding data values.

        Returns:
            dict: A dictionary where each key is a channel name, and the corresponding
                  value is a list of data points for that channel.

        Raises:
            None

        Example Usage:
            data = [
                ["t", "EAG1", "EAG4", "EAG3", "EAG2","Solenoid"],
                ["0", "0.1", "0.2", "0.3"],
                ["1", "0.4", "0.5", "0.6"],
                ["2", "0.7", "0.8", "0.9"]
            ]
            eag_dict = get_EAGs(data)
            print(eag_dict)
            # Output: {'Channel1': [0.1, 0.4, 0.7], 'Channel2': [0.2, 0.5, 0.8], 'Channel3': [0.3, 0.6, 0.9]}

        """
    chList = data[0][:]

    # Find the indices of the channels in the first row
    chIndices = [i for i, val in enumerate(data[0]) if val in chList]

    # Create an empty dictionary to store the channels and their data
    channels_data = {}

    # Loop through each channel index
    for chIndex in chIndices:
        # Get the channel name using the index
        ch = chList[chIndex]

        # Extract the data for the current channel
        channel_data = [float(data[x][chIndex]) for x in range(1, len(data))]
        # Calculate the baseline by taking the average of rows 50 to 250
        baseline = sum(channel_data[49:250]) / (250 - 50 + 1)

        # Subtract the baseline from the entire column
        baseline_subtracted_data = [value - baseline for value in channel_data]

        # Store the channel and its data in the dictionary
        channels_data[ch] = baseline_subtracted_data

    return channels_data

#'EAG2' corresponds to top left location in the chip so this is the default recording channel. others can be added.
#this will return a dictionary for all recording channels specified
def Extract_mEAG(FILE, record_channels=['EAG2'], lc=.1, hc=4, BF=True):
    # first we load the data into our variable abf
    # open a csv file containing a EAG wave
    #with open(FILE, newline='') as f:
    with open(FILE, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)

    EAGs_dict = get_EAGs(data)
    solenoid = EAGs_dict['Solenoid']
    #the stimulus data is stored in "solenoid". we need to identify when it is activated/inactivated
    sol = find_sol(solenoid)
    ni = len(sol)

    # extract a 5 second window centered on the solenoid for each wave (this means three waves per channel
    # Each key in the dictionary is either the time 't', or a channel 'EAG1', 'EAG4','EAG3','EAG2','Solenoid'
    # the values for each channel with be three lists 5 seconds in length (500 points)
    for key, value in EAGs_dict.items():
        intervals = []
        # store the channel to be processed in the variable temp
        TEMP = EAGs_dict[key]
        for i in range(0, ni):
            interval = TEMP[sol[i][0] - 50: sol[i][1] + 400]
            intervals.append(interval)
        EAGs_dict[key] = intervals

    #Apply the butterworth filter to the values of the dictionary
    if BF == True:
        for key, value in EAGs_dict.items():
            filtered_data = []
            TEMP = EAGs_dict[key]
            for i in range(3):
                fdata = butter_bandpass_filter(TEMP[i], lowcut=lc, highcut=hc, fs=1000.0, order=1)
                filtered_data.append(fdata)
            EAGs_dict[key] = filtered_data

    filtered_dict = {key: value for key, value in EAGs_dict.items() if key in record_channels}
    return filtered_dict



eag_dict= Extract_mEAG('Data/Raw_Data/062123/062123M2fA1/062123M2fA1_50ul_ylangylang1.csv', lc=.1,hc=4,BF=True)
for x in range(3):
    plt.plot(eag_dict['EAG2'][x])
    plt.show()


print(eag_dict.keys())