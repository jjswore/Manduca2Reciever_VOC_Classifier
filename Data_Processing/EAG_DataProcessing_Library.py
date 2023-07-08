import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import csv
import glob
from scipy.signal import butter, lfilter
from scipy import optimize

def find_sol(data):
    """
    This function takes an array of data, calculates the difference between consecutive 
    elements, and then finds the first indices of series where the difference is either 
    less than -0.6 or greater than 0.6.

    Arguments:
    data -- a list or array-like object of numerical values

    Returns:
    sol -- a list of tuples, where the first element of each tuple is an index where 
           the difference is greater than 0.6, and the second element is an index where 
           the difference is less than -0.6.
    """
    # Calculate the difference between consecutive elements in the array
    SolFall = np.diff(data)

    # Initialize the lists to store the indices and the flags
    NI, PI = [], []
    flag_ni, flag_pi = False, False

    # Iterate through the SolFall array starting from the 6th element
    for i, v in enumerate(SolFall[5:]):
        # If the difference is less than -0.6 and the last data point did not meet this condition
        if v < -0.6 and not flag_ni:
            # Add the current index to the NI list
            NI.append(i)
            # Set the flag to True to indicate that the current data point meets the condition
            flag_ni = True
        elif v >= -0.6:
            # If the difference is not less than -0.6, reset the flag to False
            flag_ni = False

        # If the difference is more than 0.6 and the last data point did not meet this condition
        if v > 0.6 and not flag_pi:
            # Add the current index to the PI list
            PI.append(i)
            # Set the flag to True to indicate that the current data point meets the condition
            flag_pi = True
        elif v <= 0.6:
            # If the difference is not more than 0.6, reset the flag to False
            flag_pi = False

    # Create a list of tuples where each tuple contains an element from PI and an element from NI
    test = zip(PI, NI)
    # Convert the zip object to a list
    sol = list(test)

    # Return the list of tuples
    return sol


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
def Extract_mEAG(FILE, record_channels=['EAG2'], lc=.1, hc=4, order=1, BF=True):
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

    # extract a 4 second window centered on the solenoid for each wave (this means three waves per channel
    # Each key in the dictionary is either the time 't', or a channel 'EAG1', 'EAG4','EAG3','EAG2','Solenoid'
    # the values for each channel with be three lists 5 seconds in length (400 points)

    for key, value in EAGs_dict.items():
        if key not in record_channels:
            continue
        #print(key)
        intervals = []
        # store the channel to be processed in the variable temp
        TEMP = EAGs_dict[key]
        for i in range(0, ni):
            interval = TEMP[sol[i][0] - 50: sol[i][1] + 400]

            intervals.append(interval)
        EAGs_dict[key] = intervals
        #print(EAGs_dict[key])

    #Apply the butterworth filter to the values of the dictionary
    if BF == True:
        for key, value in EAGs_dict.items():
            if key not in record_channels:
                continue
            filtered_data = []
            TEMP = EAGs_dict[key]

            for i in range(3):
                fdata = butter_bandpass_filter(TEMP[i], lowcut=lc, highcut=hc, fs=100.0, order=order)

                # Find the index of the peak (maximum absolute value)
                #print(f'this is the filtered data {fdata}')
                peak_index = np.argmax(np.abs(fdata))
                # If the peak is negative, invert the sign of all values in this trial
                if fdata[peak_index] < 0:

                    fdata = [-v for v in fdata]

                filtered_data.append(fdata)

            EAGs_dict[key] = filtered_data
    filtered_dict = {key: value for key, value in EAGs_dict.items() if key in record_channels}
    return filtered_dict

def name_con(f):
    # This will splits the basename of the file on "_" to find the concentration in the file name
    #'dateMoth#SexAntenna#_line_deliverymethod_odor_trial' 062623M1fA1_OR6KO_p_linalool_1
    tn = os.path.basename(f).split("_", 3)
    name = f'{tn[0]}{tn[3]}'
    return name

def butter_bandpass(lowcut, highcut, fs, order=1):
    # creates a butterworth bandpass filter
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')

    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=1):
    # applies a butterworth bandpass filter to some data
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def ConcCh(data):
    # concatenates the data from two channels into a single array

    waves = len(data[0])
    Ndata = np.zeros((waves, 18000))
    for ch in range(len(data)):
        for w in range(0, waves):
            try:
                Ndata[w, :] = np.concatenate((data[0][w], data[ch + 1][w]))
            except:
                pass
    return Ndata

def SaveData(data, directory, name):  # data is an numpy array
    # saves nested array of multichannel data to a

    Dir = directory
    if not os.path.isdir(Dir):
        os.makedirs(Dir)

    for key in data:
        for wave_idx, wave_data in enumerate(data[key]):
            with open(f'{Dir}/{name}_{key}_wave{wave_idx}.csv', 'w', newline='') as f:
                write = csv.writer(f)
                write.writerow(wave_data)


def namer(f):
    #removes the "." from the file name
    n = os.path.basename(f)
    tn = n.split(".")

    return tn[0]


def findCTRL(file1, folder):
    #used the find the control file.
    result = 1000000000
    ctrl = None
    for x in folder:
        #get the time difference between experiment and the control. repeat for all files in folder
        tt = abs(os.path.getmtime(file1) - os.path.getmtime(x))
        if tt < result:
            #if the difference is smaller than the result variable then replace "result" with new dif
            result = tt
            ctrl = x
            # print(ctrl)
    return ctrl

def open_wave(FILE):
    #open a csv file containing a EAG wave
    with open(FILE, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)
    f.close()
    l = data[0]
    l = list(map(float, l))
    return l

def csv_plot(FILE, NAME):
    #plot a csv file
    t = open_wave(FILE)
    # n=NAME.split("\\")
    plt.title(label=NAME, size=10)
    plt.plot(t)
    plt.show()
    plt.close()
def MinMax_Norm(xdata):
    # Normalize the data to the greatest responding odorant (Currently YlangYlang)
    # Ylangylang is used as the "control"
    mx = max(xdata)
    mn = min(xdata)
    xScaled = [((n - mn) / (mx - mn)) for n in xdata]

    # for n in xdata:
    # x=(n-mn)/(mx-mn)
    # xScaled.append(x)
    return xScaled


def Min_Normalization(minCTRL_xdata, xdata):
    """Normalize the data to the greatest responding odorant (Currently YlangYlang)
    Ylangylang is used as the "control"""

    xScaled = []
    mx = max(minCTRL_xdata)
    mn = min(minCTRL_xdata)
    for n in xdata:
        x = (n - mn) / (mx - mn)
        xScaled.append(x)
    return xScaled

def log_transform(file):
    """
    Apply log2 transformation to a numpy array.

    Parameters:
    file (ndarray): A numpy array to be transformed.

    Returns:
    ndarray: A log2 transformed numpy array.
    """
    # Get the minimum value in the array and add a small value to prevent taking the log of zero
    a = np.min(file)
    # Apply the log2 transformation to each element in the array
    log2T = [np.log2(x + a) for x in file]
    return log2T


def find_MaxIntenseWave(data):
    """
    This function takes a list of waves as input and returns the wave with the highest maximum absolute value.

    Args:
    - data: A list of lists or a list of 1D numpy arrays, each representing a wave.

    Returns:
    - The wave with the highest maximum absolute value as a list or a 1D numpy array (depending on the input type).

    Example:
    >>> data = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    >>> find_MaxIntenseWave(data)
    [7, 8, 9]
    """
    # Convert data to numpy arrays if not already
    data = [np.array(d) for d in data]

    # Find the index of the wave with the maximum absolute value
    best_index = max(range(len(data)), key=lambda i: np.abs(data[i]).max())

    # Return the wave with the highest maximum absolute value
    return data[best_index].tolist()  # Convert back to list if original data was a list


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

def Mean_Smoothing(data, window, normalized=False):
    if normalized==True:
        filtsig = [1 for x in range(len(data))]
    elif normalized ==False:
        filtsig = [0 for x in range(len(data))]
    for i in range (window+1, len(data)-window-1):
        filtsig[i] = np.mean(data[i-window:i+window])
    #windowsize = 1000*(window*2+1) / 1000
    return np.array(filtsig)

def get_subdirectories(directory):
    subdirectories = []
    for name in os.listdir(directory):
        path = os.path.join(directory, name)
        if os.path.isdir(path):
            subdirectories.append(path)
    return subdirectories

def PSD_analysis(data):
    dt = 1
    t = np.arange(9000, 18000, dt)
    n = len(t)
    fhat = np.fft.fft(data.T, n)  # Compute the FFT
    PSD = np.abs(fhat) ** 2 / (dt * n)  # Power spectrum (power per freq)
    freq = np.fft.fftfreq(n, dt)[:n//2]  # Create x-axis of frequencies in Hz
    return PSD[:n//2], freq

def FFT_analysis(data, th):
    dt = 1
    t = np.arange(0, 9000, dt)
    n = len(t)
    fhat = np.fft.fft(data.T, n)
    PSD = np.abs(fhat) ** 2 / n
    freq = ((dt * n) / 9) * np.arange(n)
    L = np.arange(1, n // 2)
    alpha, pcov = optimize.curve_fit(lambda x, a, b: a * x + b, freq[1:9], np.log(PSD[1:9]))
    if pcov[1][1] > th:
        perr = np.sqrt(np.diag(pcov))
        print('perr:', perr)
    PSD[PSD <= 2] = 0
    return pcov[1][1]

def FFT_LSTSQ_QC(df,t):
    if len(df) > 10000:
        CH1=df.T.iloc[:9000,]
        CH2=df.T.iloc[9000:-3,]
        E_List=list(df.T.columns)
        good1=[CH1.columns[x] for x in range(len(CH1.columns)) if
            (FFT_analysis(CH1[E_List[x]],th=t)) < t]
        good2=[CH2.columns[x] for x in range(len(CH2.columns)) if
            (FFT_analysis(CH2[E_List[x]],th=t)) < t]
        good=list(np.intersect1d(good1,good2))
        QCdf=df.T[good]
    else:
        CH1 = df.T.iloc[:9000, ]
        E_List = list(df.T.columns)
        good = [CH1.columns[x] for x in range(len(CH1.columns)) if
                 (FFT_analysis(CH1[E_List[x]], th=t)) < t]
        QCdf = df.T[good]
    print(len(df))
    print(len(QCdf.T))
    return(QCdf)

def find_sol1(data):#data = abf.sweepY of second channel out of three channels "channel = 1" from abf file
    SolFall = np.diff(data)
    NI = [i*.001 for i,v in enumerate(SolFall) if v < -0.3]
    PI = [i*.001 for i,v in enumerate(SolFall) if v > 0.3]
    test=zip(PI,NI)
    sol=list(test)
    return(sol)

def intensity_alignment(df):
    IDmins=list(df.iloc[:,:2600].idxmin(axis=1))
    min_idx=[int(x) for x in IDmins]
    maxSHIFT=max(min_idx)
    df_data = df.iloc[:, :-3]
    index_map = {name: i for i, name in enumerate(df_data.index)}
    shifted_df_data = df_data.apply(lambda row: row.shift(maxSHIFT - min_idx[index_map[row.name]]), axis=1)
    shifted_df_data.fillna(0, inplace=True)
    shifted_df_data = pd.concat([shifted_df_data,df.iloc[:,-3:]],axis=1)
    return shifted_df_data

def process_data(DL=[], record_channels=['EAG1','EAG2','EAG3','EAG4'], savedir=None, SUM=False, Norm='YY', Smoothen=False, LOG=False,
                 Butter=[.1, 4, 1], B_filt=True, RETURN='Save'):
    """
        Process data from DList and save the processed data to savedir.

        Parameters:
        DList (list): A list of directories containing the data to process.
        norm (str): Normalization method to use. Possible values: 'YY', 'Yes', False. Default is 'YY'. This normalizes
                    data to the strongest odorant ylangylang.
        sub (bool): Whether to subtract channel 1 from channel 2. Default is False.
        savedir (str): Directory to save processed data to. Default is ''.

        Returns: None
        """


    print('Starting Data_Processing')
    B_filt=B_filt
    Sum = SUM
    SAVEDIR = savedir
    DList=[(subdir +'/') for directory in DL for subdir in get_subdirectories(directory)]
    Normalize = Norm
    print(f'this is the DLIST: {DList}')
    for D in DList:
        print('beginning ', D)
        f1 = [f.path for f in os.scandir(D)
              if 'DS_Store' not in os.path.basename(f)]
        # seperate the data into experimental and control lists
        ctrl = [x for x in f1 if 'mineraloil' in os.path.basename(x.lower()) or 'compressedair' in os.path.basename(x.lower())]
        print(f'this is the control{ctrl}')
        exp = [x for x in f1 if 'mineraloil' not in os.path.basename(x.lower()) or 'compressedair' not in os.path.basename(x.lower())]
        YY = [x for x in f1 if 'ylangylang' in os.path.basename(x.lower())]


        # Extract the each individual wave and subtract the miniral oil control
        for data in exp:

            control = findCTRL(data, ctrl)
            RefWave = findCTRL(data, YY)
            # print(data,control)
            n = os.path.basename(data)
            print(n, 'is an experiment')
            VOC = n.split("_")[4]
            if n.split("_")[2] == 's':
                delivery = 'Syringe'
            elif n.split("_")[2] == 'p':
                delivery = 'Pipette'
            #create the directory where waves will be saved
            DIR = f'{SAVEDIR}{delivery}/{VOC}/'
            n = namer(data)
            Odor = Extract_mEAG(data,record_channels, Butter[0], Butter[1], BF=B_filt)
            CTRL = Extract_mEAG(control,record_channels, Butter[0], Butter[1],BF=B_filt)
            RW = Extract_mEAG(RefWave,record_channels, Butter[0], Butter[1],BF=B_filt)

            for x in range(0, 3):
                # subtract the control

                for key in Odor.keys():
                    Odor[key][x] = [a - b for a, b in zip(Odor[key][x], [np.mean(values) for values in zip(*[wave[:500] for wave in CTRL[key]])])]

            if Normalize == 'Yes':
                for key in Odor:
                    for x in range(3):
                        Odor[key][x] = MinMax_Norm(Odor[key][x])

            elif Normalize == 'YY':
                for key in Odor.keys():
                    RW1 = find_MaxIntenseWave(RW[key])
                    for x in range(3):
                        # normalize to ylangylang
                        Odor[key][x] = Min_Normalization(RW1, Odor[key][x])
                        # baseline the result to 0
                        avg = np.mean(Odor[key][x][0:500])
                        Odor[key][x] = [n - avg for n in Odor[key][x]]

                        if Smoothen:
                                Odor[key][x] = Mean_Smoothing(data=Odor[key][x], window=125, normalized=True)

            elif Normalize == False:
                for key in Odor:
                    for x in range(3):
                        avg = np.mean(Odor[key][x][0:500])
                        Odor[key][x] = [n - avg for n in Odor[key][x]]

            if Smoothen:
                for key in Odor:
                    for x in range(3):
                        Odor[key][x] = Mean_Smoothing(data=Odor[key][x], window=125, normalized=False)

            if LOG:
                for key in Odor:
                    for x in range(3):
                        Odor[key][x] = log_transform(Odor[key][x])

            if Sum == True:
                summed_lists = {'Summed': []}

                for i in range(3):
                    summed_wave = list(np.sum([Odor[key][i] for key in Odor.keys() if len(Odor[key]) > i], axis=0))
                    summed_lists['Summed'].append(summed_wave)

                Odor = summed_lists

            if RETURN == 'SAVE':
                SaveData(Odor, directory=DIR, name=n)


            elif RETURN == 'PLOT':
                for key in Odor.keys():
                    for w in range(3):
                        i = input('do you want to plot')
                        if i.lower() == 'yes':
                            plt.plot(Odor[key][w])
                            plt.title(f'{VOC}{key} wave{w}')
                            plt.ylim(-1.5,1.5)
                            plt.show()
                        else:
                            break
            print('finished')
