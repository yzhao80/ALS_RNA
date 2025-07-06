from numpy import loadtxt
from numpy import sort
from xgboost import XGBClassifier,XGBRFClassifier
from sklearn.model_selection import train_test_split, cross_validate, KFold
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, precision_score, recall_score
from sklearn.feature_selection import SelectFromModel
import pandas as pd
import sys
import datetime

import pyreadr

import joblib

#get the input from arg
version = sys.argv[1]
sub_dir = sys.argv[2]

filename = f"Rdata/external_training_test_{version}.rda"
current_date_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#print the process time and filename
print(f"process time: {current_date_time},filename: {filename}")

training_test = pyreadr.read_r(filename)
train1 = training_test["um_data"]
test1 = training_test["external_data"]
#let train1 feature name assign to test1
# test1.columns = train1.columns
print(f"train1 shape: {train1.shape}, test1 shape: {test1.shape}")

X1_train=train1.drop(["Row.names","ALS_status"],axis=1)
y1_train=train1["ALS_status"]
y1_train=y1_train.replace({'case':1,'control':0})

X1_test=test1.drop(["Row.names","ALS_status"],axis=1)
y1_test=test1["ALS_status"]
y1_test=y1_test.replace({'case':1,'control':0})

print(f"X1_train shape: {X1_train.shape}, X1_test shape: {X1_test.shape}")

def cross_validate_xgboost_model( X_train, y_train, n_splits=10, suppress_output=False):
    params = {"objective": "binary:logistic", 
                "eval_metric": "auc", 
                "eta":0.1, 
                "max_depth":20,
                "lambda": 0.0003, "alpha": 0.0003, "nthread" :30}
    model = XGBClassifier(**params) 


    kf = KFold(n_splits=n_splits, random_state=6, shuffle=True)
    
    #store the accuracy and auc for each fold
    accuracy_list = []
    roc_auc_list = []
    best_model = None

    for fold, (train_index, valid_index) in enumerate(kf.split(X_train)):
        X_train_fold, X_valid_fold = X_train.iloc[train_index], X_train.iloc[valid_index]
        y_train_fold, y_valid_fold = y_train.iloc[train_index], y_train.iloc[valid_index]
        model.fit(X_train_fold, y_train_fold.values)
        y_pred = model.predict(X_valid_fold)
        y_pred_proba = model.predict_proba(X_valid_fold)[:,1]
        accuracy = accuracy_score(y_valid_fold, y_pred)
        roc_auc = roc_auc_score(y_valid_fold, y_pred_proba)
        f1 = f1_score(y_valid_fold, y_pred)
        precision = precision_score(y_valid_fold, y_pred)
        recall = recall_score(y_valid_fold, y_pred)


        if roc_auc > max(roc_auc_list, default=0):
            best_model = model
        accuracy_list.append(accuracy)
        roc_auc_list.append(roc_auc)
        
        if not suppress_output:
            print(f"Fold: {fold}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, F1: {f1}, Precision: {precision}, Recall: {recall}")
        
    return {"model": best_model, "accuracy": accuracy_list, "roc_auc": roc_auc_list}

whole_results = cross_validate_xgboost_model(X1_train, y1_train, n_splits=10)
best_model_whole = whole_results["model"]
#print the average accuracy and roc_auc
print(f"Average accuracy: {sum(whole_results['accuracy'])/len(whole_results['accuracy'])}")

y1_pred = best_model_whole.predict(X1_test)
y1_pred_proba = best_model_whole.predict_proba(X1_test)[:,1]
accuracy = accuracy_score(y1_test, y1_pred)
roc_auc = roc_auc_score(y1_test, y1_pred_proba)
f1 = f1_score(y1_test, y1_pred)
precision = precision_score(y1_test, y1_pred)
recall = recall_score(y1_test, y1_pred)

print(f"Accuracy: {accuracy}, ROC AUC: {roc_auc}, F1: {f1}, Precision: {precision}, Recall: {recall}")

thresholds = sort(best_model_whole.feature_importances_)
thresholds_fit = thresholds[thresholds>0]
print(f"thresholds length = {len(thresholds_fit)}")
best_model=None
best_AUC = 0

number_of_features = []
AUCs = []
derivatives = []

best_model=None
best_AUC = 0
best_n_features = 0

best_0_50 = {"AUC": 0, "n_features": 0}
best_100_200 = {"AUC": 0, "n_features": 0}

number_of_features = []
AUCs = []

for thresh in thresholds_fit:
    
    # select features using threshold
    selection = SelectFromModel(estimator=best_model_whole, threshold=thresh, prefit=True)
    col_index = X1_train.columns[selection.get_support()]

    select_X_train = selection.transform(X1_train)
    select_X_train = pd.DataFrame(select_X_train, columns=col_index)

    print(f"thresh: {thresh}, n={select_X_train.shape[1]}")
    # train model
    results = cross_validate_xgboost_model(select_X_train, y1_train, n_splits=10, suppress_output=True)

    #get the average of accuracy and auc
    accuracy = sum(results['accuracy'])/len(results['accuracy'])
    auc = sum(results['roc_auc'])/len(results['roc_auc'])

    print(f"accuracy: {accuracy}, auc: {auc}")
    print("--------------------------------------------------")
    number_of_features.append(select_X_train.shape[1])
    AUCs.append(auc)


print(f"Best AUC: {best_AUC} \t Number of features: {best_n_features}")

#save number of features, auc in dataframe
df = pd.DataFrame({"number_of_features": number_of_features, "AUCs": AUCs})
df.to_csv(f"feature_selection_params{version}.csv")

#plot the number of features vs auc
import matplotlib.pyplot as plt
plt.plot(number_of_features, AUCs)
#add the 0-50 and 100-200 best Auc vertical line
plt.xlabel('Number of features')
plt.ylabel('AUC')
plt.savefig(f'feature_selection_{version}.png')



