import random
import pandas as pd
import numpy as np
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingGridSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import make_scorer
from sklearn.metrics import auc

import pickle
from sklearn.metrics import recall_score

def name_con(f):
    n=os.path.basename(f)
    tn=n.split("_")
    if tn[1].lower()=='100':
        x='1:100'
    elif tn[1].lower()=='1k':
        x='1:1k'
    elif tn[1].lower()=='10k':
        x='1:10k'
    else:
        x=''
    name=x+' '+tn[2]
    return name

def AUC(DF, WF,xPeak1,xPeak2):
    AUC= auc(abs(DF[WF][xPeak1:xPeak2].index), abs(DF[WF][xPeak1:xPeak2]))
    return AUC

def PCA_Constructor(data, PC=3):
    wave_df=data.iloc[:,:-3].convert_dtypes(float).astype(float)
    scaled_data=pd.DataFrame(StandardScaler().fit_transform(wave_df),columns = wave_df.columns,
                         index=wave_df.index)
    PCA_set=PCA(n_components=PC)
    PCA_Results=PCA_set.fit_transform(scaled_data)
    PCA_DF = pd.DataFrame(data = PCA_Results, index=data.index)
    for pc in range(PC):
        PCA_DF.rename({pc:f'PC {(pc+1)}'}, axis=1, inplace=True)
    return PCA_DF, PCA_set, scaled_data

def TT_Split(DF, t=.7): #splits into train and test data by antenna so that a given percentage of antennas are used for training
    List=[x for x in DF['date'].unique()]
    Shuffled=random.sample(List, k=len(List))

    TrainDF = DF.loc[DF['date'].isin(Shuffled[0:round(len(Shuffled)*t)])]
    TestDF = DF.loc[DF['date'].isin(Shuffled[round(len(Shuffled)*t):])]
    TrainL = TrainDF['label']
    TestL = TestDF['label']
    return(TrainDF.drop(['label','date','concentration'], axis=1),
            TestDF.drop(['label','date','concentration'], axis=1),
           TrainL, TestL)

def RFC_GridSearch(data, concentration, odors):
    Analysis_data=pd.concat(data)

    data_df=Analysis_data[Analysis_data['concentration'].str.contains(concentration)]
    data_df=data_df[data_df['label'].str.contains(odors)]

    print('splitting data')
    train_features, test_features, train_labels, test_labels = TT_Split(data_df, .7)

    n_estimators = [10]
    max_features = ["sqrt", "log2"]
    max_depth = list(np.arange(10, 120, step=10))
    max_leaf_nodes = list(np.arange(10,510, step=50))
    min_samples_split = [2,4,6,8]
    min_samples_leaf = [1,2,3,4,8]
    max_samples=[.1,.2,.3,.4,.5,.6,.7,.8]
    bootstrap = [True]
    ###### GRID SEARCH ########

    param_grid = {
        "n_estimators":n_estimators,
        "max_features": max_features,
        "max_depth": max_depth,
        "min_samples_split": min_samples_split,
        'max_leaf_nodes':max_leaf_nodes,
        "min_samples_leaf": min_samples_leaf,
        "max_samples": max_samples,
        "bootstrap":bootstrap}



    acc_scorer = make_scorer(accuracy_score)
    print("beginning grid search")
    GRID_cv = HalvingGridSearchCV(RandomForestClassifier(),  param_grid, scoring=acc_scorer, cv=15, n_jobs=7, min_resources=14, error_score='raise')
    GRID_cv.fit(train_features, train_labels)
    gbp=GRID_cv.best_params_
    gbs=GRID_cv.best_score_

    rfc = RandomForestClassifier(n_estimators=1000, min_samples_split=gbp['min_samples_split'],
                                 min_samples_leaf=gbp['min_samples_leaf'],max_features=gbp['max_features'],
                                 max_depth=gbp['max_depth'],bootstrap=gbp['bootstrap'],
                                 max_leaf_nodes=gbp['max_leaf_nodes'],
                                 max_samples=gbp['max_samples'], oob_score=True, n_jobs=7)

    print('Classifier Parameters Found', ": ", rfc)
    return rfc, gbp, gbs

def RF_Testing(data, concentration, odors, classifier, P):
    """
    Trains a Random Forest classifier on the input data and returns accuracy score, confusion matrix,
    predictions, probabilities, predicted classes and log probabilities.

    Args:
    - data (pandas.DataFrame): Input data
    - concentration (str): Concentration to use for analysis
    - odors (str): Odor label to use for analysis
    - classifier (sklearn.ensemble.RandomForestClassifier): Random Forest classifier object

    Returns:
    - results_dict (dict): A dictionary containing the following keys:
      - 'accuracy_score': Accuracy score
      - 'confusion_matrix': Confusion matrix
      - 'predictions': Predictions
      - 'probabilities': Probabilities
      - 'predicted_classes': Predicted classes
      - 'log_probabilities': Log probabilities
    """

    # Concatenate the input data
    Analysis_data=pd.concat(data)

    # Get the data for the given concentration and odor label
    title=concentration
    data_df=Analysis_data[Analysis_data['concentration'].str.contains(title)]
    data_df=data_df[data_df['label'].str.contains(odors)]

    # Split the data into training and testing sets
    print('splitting data')
    train_features, test_features, train_labels, test_labels = TT_Split(data_df, .7)
    params = classifier.get_params()
    params['n_estimators']=1000
    classifier=RandomForestClassifier(**params)
    print(classifier.get_params(deep=True))
    # Train the model on training data
    classifier.fit(train_features, train_labels)

    # Get predictions and probabilities on the testing set
    predictions = classifier.predict(test_features)
    probabilities = classifier.predict_proba(test_features)
    logprob= classifier.predict_log_proba(test_features)

    # Get the confusion matrix and accuracy score
    CM=confusion_matrix(test_labels,predictions, labels=classifier.classes_).astype(float)
    acc_score=balanced_accuracy_score(predictions, test_labels)
    r_score = recall_score(test_labels, predictions, zero_division=0, average='macro', labels=[P])
    print('Accuracy: ', acc_score, '\n')

    # Create a dictionary of results
    results_dict = {'classifier':classifier,
                    'accuracy_score': acc_score,
                    'sensitivity': r_score,
                    'confusion_matrix': CM,
                    'true classes':test_labels,
                    'predictions': predictions,
                    'probabilities': probabilities,
                    'predicted_classes': classifier.classes_,
                    'log_probabilities': logprob}
    return results_dict

def SVM_GridSearch(data, concentration, odors):
    """
    Perform a grid search to optimize SVM hyperparameters.

    Args:
    - data (List[pd.DataFrame]): A list of pandas dataframes, each containing the data to be analyzed
    - concentration (str): The concentration of the odor stimuli to be analyzed
    - odors (str): The label of the odor stimuli to be analyzed
    - P (str): The positive class label for computing recall score

    Returns:
    - clf (svm.SVC): The optimized SVM classifier
    - gbp (Dict[str, Any]): The best set of hyperparameters found by grid search
    - gbs (float): The best score found by grid search
    """

    # Concatenate data into a single dataframe
    Analysis_data = pd.concat(data, axis=1)

    # Filter data based on concentration and odor label
    data_df = Analysis_data[(Analysis_data['concentration'].str.contains(concentration)) &
                            (Analysis_data['label'].str.contains(odors))]

    # Split data into train and test sets
    print('Splitting data...')
    train_features, test_features, train_labels, test_labels = TT_Split(data_df, .7)

    # Set hyperparameters to search over
    kernel = ['rbf']
    C = [1, 2, 2, 5, 3, 4, 5, 7, 7.5, 8]
    degree = [0, 0.01, 0.05, 0.1, 0.5]
    gamma = ['scale', 'auto', 0.1, 0.2, 0.5]
    coef0 = [0, 0.05, 0.1, 0.2]

    # Create parameter grid
    param_grid = {
        "kernel": kernel,
        "C": C,
        "degree": degree,
        "gamma": gamma,
        "coef0": coef0
    }

    print("Beginning grid search...")
    # Perform grid search
    GRID_cv = GridSearchCV(
        svm.SVC(),
        param_grid,
        scoring='accuracy',
        n_jobs=-1,
        error_score='raise',
        cv=15,
        verbose=1
    )
    GRID_cv.fit(train_features, train_labels)

    # Extract best hyperparameters and score
    gbp = GRID_cv.best_params_
    gbs = GRID_cv.best_score_

    # Create optimized classifier
    clf = svm.SVC(
        kernel=gbp['kernel'],
        C=gbp['C'],
        degree=gbp['degree'],
        gamma=gbp['gamma'],
        coef0=gbp['coef0']
    )

    print(f"Best parameters found: {gbp}")
    print(f"Best score found: {gbs}")
    print(f"Optimized classifier: {clf}")

    return clf, gbp, gbs

def SVM_Testing(data, concentration, odors, classifier, P):
    """
    Applies SVM classification on the provided data.

    Args:
    - data: pandas.DataFrame, input dataset
    - concentration: str, the concentration of odorant to classify
    - odors: str, the odorants to classify
    - classifier: sklearn.svm.SVC, the classifier to use
    - P: str, the positive class

    Returns:
    - results_dict: dict, a dictionary containing the results of the classification
    """

    # Filter the dataset based on concentration and odorants
    # Concatenate the input data
    Analysis_data=pd.concat(data)

    # Get the data for the given concentration and odor label
    title=concentration
    data_df=Analysis_data[Analysis_data['concentration'].str.contains(title)]
    data_df=data_df[data_df['label'].str.contains(odors)]

    # Split the dataset into training and testing sets
    print('Splitting data...')
    train_features, test_features, train_labels, test_labels = TT_Split(data_df, .7)

    # Train the classifier
    print('Training model...')
    classifier.fit(train_features, train_labels)

    # Predict the response for test dataset
    print('Predicting classes...')
    predictions = classifier.predict(test_features)

    # Calculate accuracy and recall
    acc_score = balanced_accuracy_score(predictions, test_labels)
    r_score = recall_score(test_labels, predictions, zero_division=0, average='macro', labels=[P])

    # Compute the confusion matrix
    CM = confusion_matrix(test_labels, predictions, labels=classifier.classes_).astype(float)

    # Compute the decision function values
    probabilities = classifier.decision_function(test_features)

    # Print the results
    print('Accuracy: ', acc_score, '\n')
    print('Sensitivity: ', r_score, '\n')

    # Create a dictionary containing the results
    results_dict = {'classifier': classifier,
                    'accuracy_score': acc_score,
                    'sensitivity': r_score,
                    'confusion_matrix': CM,
                    'true classes': test_labels,
                    'predictions': predictions,
                    'probabilities': probabilities,
                    'predicted_classes': classifier.classes_}

    return results_dict

def LogR_GridSearch(data, concentration, odors):
    Analysis_data = pd.concat(data, axis=1)

    data_df = Analysis_data[Analysis_data['concentration'].str.contains(concentration)]
    data_df = data_df[data_df['label'].str.contains(odors)]

    print('splitting data')
    train_features, test_features, train_labels, test_labels = TT_Split(data_df, .7)
    #create an itereable for cross validation in hyperparameter tuning...
    #folds[[x for TT_Split(train_features, .7) in range(5)]]

    C = [.1, 1, 2, 2, 5, 3, 4, 5, 7, 7.5, 8, 9, 10]
    fit_intercept = [True, False]
    solver = ['liblinear']
    penalty = ['l2']
    l1_ratio = [.1, .2, .3, .4, .5, .6, .7, .8, .9]

    param_grid = {
        "C": C,
        "fit_intercept": fit_intercept,
        "solver": solver,
        "penalty": penalty}
    # 'l1_ratio':l1_ratio}

    GRID_cv = GridSearchCV(LogisticRegression(), param_grid, cv=15, scoring='accuracy', n_jobs=7,
                           error_score='raise')  # , verbose=3)
    # print(train_labels)
    GRID_cv.fit(train_features, train_labels)
    print('TRAIN SCORE: ', GRID_cv.score(train_features, train_labels))
    gbp = GRID_cv.best_params_
    gbs = GRID_cv.best_score_

    print("Best score:", gbs)
    print("Best params:", gbp)

    LR = LogisticRegression(solver=gbp['solver'], C=gbp['C'],
                            penalty=gbp['penalty'], fit_intercept=gbp['fit_intercept'])
    return LR, gbp, gbp


def LogR_Testing(data, clf, odors, concentration, P): #F are the features you want to use for SVM.
    Analysis_data=pd.concat(data)

    # Get the data for the given concentration and odor label
    data_df=Analysis_data[Analysis_data['concentration'].str.contains(concentration)]
    data_df=data_df[data_df['label'].str.contains(odors)]

    train_features, test_features, train_labels, test_labels = TT_Split(data_df, .7)
    #Create a svm clf # Linear Kernel

    #Train the model using the training sets
    print('training model...')
    clf.fit(train_features,train_labels)


    #Predict the response for test dataset
    print('predicting classess...')
    predictions = clf.predict(test_features)

    #find accuracy of the classifier
    acc_score=balanced_accuracy_score(predictions, test_labels)
    r_score=recall_score(test_labels, predictions, zero_division=0, average='macro',labels=[P])
    #print(predictions)
    #compute the confusion matrix
    CM=confusion_matrix(test_labels,predictions, labels=clf.classes_).astype(float)
    probabilities=clf.predict_proba(test_features)
    print('Accuracy: ', acc_score, '\n')
    print('Sensitivity: ', r_score, '\n')

    results_dict = {'classifier':clf,
                    'accuracy_score': acc_score,
                    'sensitivity':r_score,
                    'confusion_matrix': CM,
                    'true classes': test_labels,
                    'predictions': predictions,
                    'probabilities': probabilities,
                    'predicted_classes': clf.classes_}
    return results_dict

def RFmodel(concentrations, data, odor, PosL, repeats):
    """
    Trains and evaluates a Random Forest classifier on the provided dataset for each concentration of a given odor.

    Args:
        concentrations (list): A list of concentrations of the odor to be analyzed.
        data (pd.DataFrame): A pandas DataFrame containing the training data.
        odor (str): The name of the odor being analyzed.

    Returns:
        list: A nested list containing the classification performance metrics for each iteration of the Random Forest classifier
            for each concentration of the given odor.
    """
    results = []
    for conc in [concentrations]:
        print(f"Beginning analysis for {odor} at {conc} concentration")
        classifier, params, best_score = RFC_GridSearch(data=data, concentration=conc, odors=odor)
        print(f"Best classifier for {odor} at {conc} concentration is {classifier}")
        print(f"Building Random Forest for {odor} at {conc} concentration")
        results.append([RF_Testing(data=data, classifier=classifier, concentration=conc, odors=odor, P=PosL)
                        for _ in range(repeats)])
        print(f"Finished analysis for {odor} at {conc} concentration")
    return results

def LRmodel(concentrations, data, odor, PosL, repeats):
    """
    Trains and evaluates a logistic regression classifier on the provided dataset for each concentration of a given odor.

    Args:
        concentrations (list): A list of concentrations of the odor to be analyzed.
        data (pd.DataFrame): A pandas DataFrame containing the training data.
        odor (str): The name of the odor being analyzed.

    Returns:
        list: A nested list containing the classification performance metrics for each iteration of the logistic regression classifier
            for each concentration of the given odor.
    """
    results = []
    for conc in [concentrations]:
        print(f"Beginning analysis for {odor} at {conc} concentration")
        classifier, params, best_score = LogR_GridSearch(data=data, concentration=conc, odors=odor)
        print(f"Best classifier for {odor} at {conc} concentration is {classifier}")
        print(f"Building Random Forest for {odor} at {conc} concentration")
        results.append([LogR_Testing(data=data, clf=classifier, concentration=conc, odors=odor, P=PosL)
                        for _ in range(repeats)])
        print(f"Finished analysis for {odor} at {conc} concentration")
    return results

def SVMmodel(concentrations, data, odor, PosL, repeats):
    """
    Trains and evaluates a Support Vector Machine classifier on the provided dataset for each concentration of a given odor.

    Args:
        concentrations (list): A list of concentrations of the odor to be analyzed.
        data (pd.DataFrame): A pandas DataFrame containing the training data.
        odor (str): The name of the odor being analyzed.

    Returns:
        list: A nested list containing the classification performance metrics for each iteration of the Support Vector Machine classifier
            for each concentration of the given odor.
    """
    results = []
    for conc in [concentrations]:
        print(f"Beginning analysis for {odor} at {conc} concentration")
        classifier, params, best_score = SVM_GridSearch(data=data, concentration=conc, odors=odor)
        print(f"Best classifier for {odor} at {conc} concentration is {classifier}")
        print(f"Building SVM for {odor} at {conc} concentration")
        results.append([SVM_Testing(data=data, classifier=classifier, concentration=conc, odors=odor, P=PosL)
                        for _ in range(repeats)])
        print(f"Finished analysis for {odor} at {conc} concentration")
    return results

def pickle_Saver(savedir,ext,data):
    if not os.path.isdir(savedir):
        os.makedirs(savedir)
    reader=open(f'{savedir}{ext}.pickle','wb')
    pickle.dump(obj=data,file=reader)
    reader.close()



