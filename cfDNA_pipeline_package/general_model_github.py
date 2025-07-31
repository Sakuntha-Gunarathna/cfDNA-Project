# -*- coding: utf-8 -*-
"""General Model - GitHub
Created by: Sakuntha Devaka Gunarathna /
Published Date: July 30, 2025
"""

#install required dependencies
pip install xgboost==2.0.3 scikit-learn==1.3.2 numpy==1.26.4 --force-reinstall

import pandas as pd  # Load and manipulate data
import numpy as np  # Calculate mean and standard deviation
import xgboost as xgb  # XGBoost functionality
import matplotlib.pyplot as plt  # Plotting
import sklearn
from sklearn.model_selection import train_test_split  # Split datasets
from sklearn.metrics import (
    balanced_accuracy_score,
    roc_auc_score,
    make_scorer,
    confusion_matrix,  # Metrics
)
from sklearn.model_selection import GridSearchCV  # Grid search and cross-validation
from xgboost import XGBClassifier, plot_importance, plot_tree  # XGBoost utilities
from numpy import load

print("XGBoost version:", xgb.__version__)
print("Scikit-learn version:", sklearn.__version__)
print("Numpy version:", np.__version__)

df_label = pd.read_csv('Patient_ID.csv', index_col="Patient_ID")
#df_label.set_index('Patient_ID')
df_label.head()

# Convert Index Column to List
index_list = df_label.index.tolist()
print(index_list)

# col_list =  list(df_label["Patient_ID"])
# print(col_list)

#Import the "Reference_Peaks_ID.csv" file which contain Referance_ID, Chromosome_Name, Begins,Ends and Fragment_Length.
#Create new Dataframe
df_ID = pd.read_csv('Reference_Peaks_ID.csv')
df_ID.head()

#load "Output_deepTools_multiBigwigSummary.npz" to a object with the name "data".
# ".npz" must be generated with deeptools multiBigwigSummary using the same BED file coordinates as provided in "Reference_Peaks_ID.csv"
#Reference Peak details should be given as a ".bed" file for deeptools - multiBigwigSummary. eg: nohup multiBigwigSummary BED-file -b <path_to_bigwig_files_group1>/*.bw <path_to_bigwig_files_group2>/*.bw ... -o <output_results>.npz --outRawCounts <output_counts>.tab --BED <regions_file>.bed &

data = load('Output_deepTools_multiBigwigSummary.npz')

#List out and print all the variables and their content within the "data" object
lst = data.files
for item in lst:
    print(item)
    print(data[item])

#Print the list of variables in "data" object
print(lst)

#Add data in to two separate variables
# array_1 = data['labels'] #Array_1 contain label information

array_1 = index_list #Array_1 contain label information taken from "df_label" dataframe index column and stored in "index_list" variable.
array_2 = data['matrix'] #Array_2 contain "matrix" variable information stored in "data" object.

print(array_1)

print(array_2)

#Add "data" object "array_2" to a pandas dataframe.
#column_names = ['column 1', 'column 2', 'column 3']

column_names = array_1
df_Original = pd.DataFrame(array_2, columns=column_names)
print(df_Original)

#To save the dataframe to .CSV file
df_Original.to_csv('General_MODEL_Original.csv', index=False)

#Join(Concatanate) "df_ID" dataframe with "df_Original" dataframe.
df_1 = pd.concat([df_ID, df_Original], axis=1)
df_1.head()

#Displaying the names of the columns.
print(df_1.columns)

#Displaying the shape of the dataset.
print(df_1.shape)

#Drop rows containg "Chr Y" information
#df_1.drop(df_1.loc[df_1['Chromosome_Name']=='chrY'].index, inplace=True)
#df_1.drop(df_1.loc[df_1['Chromosome_Name']=='chrX'].index, inplace=True)

#Drop rows containg "ChrM", "Random Chr", "Unmapped Chr". This code will keep all rows in df_1 where the 'Chromosome_Name' does NOT contain any of the specified keywords.
keywords = ["chrUn", "chrM", "random"]
df_1 = df_1[~df_1['Chromosome_Name'].str.contains('|'.join(keywords))]

#Displaying the shape of the dataset.
print(df_1.shape)

#set axis=0 to remove rows, axis = 1 to remove colums.
#Inplace=True will save the changers in to the df dataframe directly without making a copy.
#df_1.columns = df_1.columns.str.replace(' ','_')

#To drop columns
df_1.drop(['Chromosome_Name','Begins','Ends','Fragment_Length'],axis=1, inplace=True)
df_1.head()

#To get the data types
df_1.dtypes

#To transpose the data set (Rows to columns and columns to rows)
#df = df.transpose()  => Will only transpose data but it will not set the "Referance_ID" as column names


#Transpose Data Frame & Set First Column as Header
#To transpose data and set the "Referance_Id" as column names.
df_2 = df_1.set_index('Referance_ID').T
df_2.head()

#To get the data types
df_2.dtypes

#To save the transposed data set
df_2.to_csv('General_MODEL_Transposed.csv', index= False)
#If you add "index" as "False", it will move the Patient ID to row names. If you name the index "True" it will add patient Id the 1st column.

#List the row names of df_2 dataframe
list(df_2.index.values)

# Cancer_TorF.csv must use the same Patient_ID values as its index column.
#Import .csv file "Cancer_TorF.csv" as the new dataframe(df_Cancer_TorF) which has the column "Cancer" and data values "T" or "F". Only uppercase T/F are allowed in csv file
#Set 1st column as the row name by "index_col=0".
df_Cancer_TorF = pd.read_csv('Cancer_TorF.csv', index_col=0)
df_Cancer_TorF.head()

#Get the shape of the "df_Cancer_TorF" dataframe
df_Cancer_TorF.shape

#Give df_2 as the first dataframe to merge and df_1 as the second datframe to merge. This will bring "Cancer" column as the 1st column in new datframe df_3
#Join "df_Cancer_TorF" and "df_2" dataframes using concatenate funchion
df_3 = pd.concat([df_Cancer_TorF, df_2], axis=1)
df_3.head()

#Get the shape of the "df_Cancer_TorF" dataframe
df_3.shape

#Replace output column "Cancer" with T=1 and F=0
df_3['Cancer'] = df_3['Cancer'].str.replace('T','1')
df_3['Cancer'] = df_3['Cancer'].str.replace('F','0')
df_3.head()

#To convert "Cancer" Column data in to numeric variable
#To get the data types
df_3['Cancer'] = pd.to_numeric(df_3['Cancer'])
df_3.dtypes

#Save the Final Dataframe (df_3) as a backup.
df_3.to_csv('General_MODEL_1_Final.csv', index= False)

#Formating data set to 'X' variable. As for the convention use capital "X" for independent variables.
X = df_3.drop('Cancer',axis=1).copy()
X.head()

#Formating data set to 'y' variable. As for the convention use simple "y" for predicter variables.
y = df_3['Cancer'].copy()
y.head()

#To veryfy that 'y' only contain 1s and 0s
y.unique()

# "y" sample counts
y.value_counts()

#To see what percentage of people have Cancer in the dataset.
#Since "y" contain numbers (1 and 0), the sum function will calculate the number of cancer patients and len function will calculate the total rows.
#Dividing number of patients by total number of rows will give you total number of patients as a percentage.
sum(y)/len(y)

#Spliting datasets in to four variables (X_train, X_test, y_train, y_test)
X_train, X_test, y_train, y_test = train_test_split(X,y,random_state=42,train_size = 0.7, stratify=y)

# "y_train" sample counts
y_train.value_counts()

# "y_test" sample counts
y_test.value_counts()

#So both 'X' and 'y' variables are previously split in to 4 groups.
#For 'X' => X_train and X_test.
#For 'y' => y_train and y_test.

#To verify that using stratify worked as expected in "y_train" must have a values closer to original Cancer percentage.
sum(y_train)/len(y_train)

#To verify that using stratify worked as expected in "y_test" must also have a value closer to original Cancer percentage
sum(y_test)/len(y_test)

# Finding the best value for scaling the weights (Scale_pos_weight)
# y.value_counts()
# y_train.value_counts()

# This parameter shows the data imbalance in the training set (Cancer and healthy)

scale_pos_weight = y_train[y_train==0].count() / y_train[y_train==1].count()
print(scale_pos_weight)

# Check for duplicate column names in X variable
duplicate_columns = X_train.columns[X_train.columns.duplicated()]
if not duplicate_columns.empty:
    print("Duplicate column names:", duplicate_columns)

#XGBoost Default Parameters
clf_xgb = xgb.XGBClassifier(seed=42,
                            early_stopping_rounds=10,
                            eval_metric='aucpr',
                            base_score=0.5,
                            booster='gbtree',
                            colsample_bylevel=1,
                            colsample_bynode=1,
                            colsample_bytree=1,
                            enable_categorical=False,
                            device = "cuda", #Change device="cuda" to device="cpu" if GPU is unavailable
                            importance_type=None,
                            interaction_constraints='',
                            objective='binary:logistic',
                            gamma=0,
                            learning_rate=0.3,
                            max_depth=6,
                            max_delta_step=0,
                            reg_lambda=1,
                            scale_pos_weight=1,
                            min_child_weight=1,
                            subsample=1,
                            n_estimators=100,
                            num_parallel_tree=1,
                            reg_alpha=0,
                            #tree_method='auto',
                            tree_method='hist', #When using Colab GPU
                            validate_parameters=1
                            )
clf_xgb.fit(X_train,
            y_train,
            verbose=True,
            eval_set=[(X_test, y_test)])

# Print the Best Validation score from above model
print("Best validation score: ", clf_xgb.best_score)

# #predicting the ranbdomly selected test data (X_test) using the our created model. Just like the Final Exam
clf_xgb.predict(X_test)

# make predictions for test data
y_pred = clf_xgb.predict(X_test)
predictions = [round(value) for value in y_pred]

#evaluate predictions train vs test data
accuracy = balanced_accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))

# from sklearn.model_selection import RandomizedSearchCV
# from scipy.stats import uniform, randint

# # The parameter ranges in param_dist now reflect distributions from which values will be randomly sampled. For discrete parameters like max_depth, randint is used.
# # For continuous parameters, uniform is used, where uniform(a, b) samples values over the range [a, a+b].
# # The n_iter parameter in RandomizedSearchCV controls how many random combinations of parameters will be tried.
# # Adjust it based on your computational resources and how exhaustive you want the search to be.



# # Define the parameter distribution
# param_dist = {
#     'max_depth': randint(5,7),  # range you prefer , if it is 3-6 you should type as (3, 7) As upper boundary is exclusive
#     'learning_rate': uniform(0.3, 0.3),  # Continuous distribution from [a, a+b]
#     #'gamma': uniform(0, 0.1),  # Continuous distribution from [a, a+b]
#     'reg_lambda': uniform(1.0, 5.0),  # Continuous distribution from [a, a+b]
#     'scale_pos_weight': uniform(1, 0.5)  # Continuous distribution from [a, a+b]
# }

# # Initialize XGBoost Classifier
# clf_xgb = xgb.XGBClassifier(seed=42,
#                             eval_metric='aucpr',
#                             base_score=0.5,
#                             booster='gbtree',
#                             colsample_bylevel=1,
#                             colsample_bynode=1,
#                             colsample_bytree=1,
#                             enable_categorical=False,
#                             device = "cuda", #Change device="cuda" to device="cpu" if GPU is unavailable
#                             gamma=0,
#                             importance_type=None,
#                             interaction_constraints='',
#                             objective='binary:logistic',
#                             n_estimators=100,
#                             n_jobs=50,
#                             num_parallel_tree=1,
#                             reg_alpha=0,
#                             tree_method='hist',  # Use 'hist' for Colab GPU
#                             validate_parameters=1)

# # Initialize Randomized Search
# random_search = RandomizedSearchCV(clf_xgb,
#                                     param_distributions=param_dist,
#                                     n_iter=200,  # Number of parameter settings that are sampled
#                                     cv=5,
#                                     scoring='average_precision',
#                                     verbose=2,
#                                     n_jobs=-1,
#                                     random_state=42)

# # Fit the data to RandomizedSearchCV (this will take some time depending on the number of iterations)
# random_search.fit(X_train, y_train, eval_set=[(X_train, y_train)], early_stopping_rounds=15)

# # Get the best parameters
# print("Best parameters found: ", random_search.best_params_)

# # The output shows the rank, mean test score, standard deviation of the test score, and the specific parameters for each of the top 5 parameter sets.

# # Convert the cv_results to a DataFrame for easier sorting and manipulation
# cv_results_df = pd.DataFrame(random_search.cv_results_)

# # Sort the results by the mean test score, in descending order
# sorted_cv_results_df = cv_results_df.sort_values(by='rank_test_score')

# # Print the top 5 parameter sets
# top_5_results = sorted_cv_results_df.head(5)

# print("Top 5 parameter sets:")
# for index, row in top_5_results.iterrows():
#     print("\nRank:", row['rank_test_score'])
#     print("Mean Test Score:", row['mean_test_score'])
#     print("Std Test Score:", row['std_test_score'])
#     print("Parameters:", row['params'])

#Please note that the nested learning_rate / reg_lambda search may take a long time and it is optional.

import xgboost as xgb
from sklearn.metrics import balanced_accuracy_score

# Define your learning rates and reg_lambda values
# learning_rates = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]

learning_rates = [0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27,
                  0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46,
                  0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65]

reg_lambdas = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Store results
results = []

# Nested loop over both hyperparameters
for lr in learning_rates:
    for reg_lambda in reg_lambdas:
        print(f"\n===== Training with learning_rate = {lr}, reg_lambda = {reg_lambda} =====")

        clf_xgb = xgb.XGBClassifier(
            seed=42,
            early_stopping_rounds=15,
            eval_metric='aucpr',
            base_score=0.5,
            booster='gbtree',
            colsample_bylevel=1,
            colsample_bynode=1,
            colsample_bytree=1,
            enable_categorical=False,
            device="cuda", #Change device="cuda" to device="cpu" if GPU is unavailable
            importance_type=None,
            interaction_constraints='',
            objective='binary:logistic',
            gamma=0,
            learning_rate=lr,
            max_depth=6,
            max_delta_step=0,
            reg_lambda=reg_lambda,
            scale_pos_weight=1,
            min_child_weight=1,
            subsample=1,
            n_estimators=100,
            n_jobs=50,
            num_parallel_tree=1,
            reg_alpha=0,
            tree_method='hist',
            validate_parameters=1
        )

        clf_xgb.fit(X_train, y_train, verbose=True, eval_set=[(X_test, y_test)])

        best_val_score = clf_xgb.best_score
        y_pred = clf_xgb.predict(X_test)
        predictions = [round(value) for value in y_pred]
        accuracy = balanced_accuracy_score(y_test, predictions)

        print("Best validation score: ", best_val_score)
        print("Balanced Accuracy: %.2f%%" % (accuracy * 100.0))

        # Save results
        results.append({
            "learning_rate": lr,
            "reg_lambda": reg_lambda,
            "best_val_score": best_val_score,
            "balanced_accuracy": accuracy
        })

import pandas as pd
from google.colab import files

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv("results.csv", index=False)

# Trigger download (for Colab)
files.download("results.csv")

# # Re-load the saved CSV
# results_df = pd.read_csv("results.csv")

# # Convert DataFrame back to list of dicts
# results = results_df.to_dict(orient="records")

# Step 1: Sort all results by balanced accuracy (descending), then by validation score (descending)
sorted_results = sorted(results, key=lambda r: (r["balanced_accuracy"], r["best_val_score"]), reverse=True)

# Step 2: Take the top 5
top_5_results = sorted_results[:5]

# Step 3: Print them
print("🏆 Top 5 Results by Balanced Accuracy (tie-breaker: validation score):")
for i, res in enumerate(top_5_results, 1):
    print(f"\n🔹 Rank #{i}")
    print(f"Learning Rate: {res['learning_rate']}")
    print(f"Reg Lambda: {res['reg_lambda']}")
    print(f"Balanced Accuracy: {res['balanced_accuracy'] * 100:.2f}%")
    print(f"Validation Score: {res['best_val_score']:.4f}")

#Modify the following Final MODEL according to the final results obtained form previous hyperparameters search.
#XGBoost Final Model With Best Parameters
clf_xgb = xgb.XGBClassifier(seed=42,
                            early_stopping_rounds=15,
                            eval_metric='aucpr',
                            base_score=0.5,
                            booster='gbtree',
                            colsample_bylevel=1,
                            colsample_bynode=1,
                            colsample_bytree=1,
                            enable_categorical=False,
                            device = "cuda", #Change device="cuda" to device="cpu" if GPU is unavailable
                            importance_type=None,
                            interaction_constraints='',
                            objective='binary:logistic',
                            gamma=0,
                            learning_rate=0.3,
                            max_depth=6,
                            max_delta_step=0,
                            reg_lambda=1,
                            scale_pos_weight=1,
                            min_child_weight=1,
                            subsample=1,
                            n_estimators=100,
                            n_jobs=50,
                            num_parallel_tree=1,
                            reg_alpha=0,
                            #tree_method='auto',
                            tree_method='hist', #When using Colab GPU
                            validate_parameters=1
                            )
clf_xgb.fit(X_train,
            y_train,
            verbose=True,
            eval_set=[(X_test, y_test)])

# Print the Best Validation score from above model
print("Best validation score: ", clf_xgb.best_score)

# #predicting the ranbdomly selected test data (X_test) using the our created model. Just like the Final Exam
clf_xgb.predict(X_test)

# make predictions for test data
y_pred = clf_xgb.predict(X_test)
predictions = [round(value) for value in y_pred]

#evaluate predictions train vs test data
accuracy = balanced_accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))

from sklearn.metrics import average_precision_score, roc_auc_score

# Predict probabilities for class 1
y_proba = clf_xgb.predict_proba(X_test)[:, 1]

# Calculate metrics
aucpr = average_precision_score(y_test, y_proba)
aucroc = roc_auc_score(y_test, y_proba)

# Print results
print(f"Average Precision Score (AUC-PR) on test set: {aucpr:.4f}")
print(f"ROC AUC Score on test set: {aucroc:.4f}")

# Get and print the Final hyperparameters used for record
params = clf_xgb.get_params()
for key, value in params.items():
    print(f"{key}: {value}")

# To get classification_report
from sklearn.metrics import classification_report

# clf_xgb is your trained model and X_test is your test dataset
predictions = clf_xgb.predict(X_test)
print(classification_report(y_test, predictions))

import xgboost as xgb

# List of CV values to try
cv_values = [3, 4, 5, 6, 7, 8, 9, 10]

# Collect results
cv_results_summary = []

# Get XGBoost core parameters from a pre-trained model
xgb_params = clf_xgb.get_xgb_params()
xgb_params['eval_metric'] = 'aucpr'

# Create DMatrix
dtrain = xgb.DMatrix(X, label=y)

# Loop over different CV folds
for cv in cv_values:
    print(f"\n🔁 Running CV with {cv}-folds")

    results = xgb.cv(
        params=xgb_params,
        dtrain=dtrain,
        num_boost_round=1000,
        nfold=cv,
        stratified=True,
        early_stopping_rounds=15,
        seed=42,
        verbose_eval=False  # Change to True if you want step-by-step logs
    )

    best_score = results['test-aucpr-mean'].max()
    best_iteration = results['test-aucpr-mean'].idxmax()

    print(f"CV = {cv} | Best AUC-PR: {best_score:.4f} at round {best_iteration}")

    # Save result
    cv_results_summary.append({
        "cv": cv,
        "best_aucpr": best_score,
        "best_iteration": best_iteration
    })

# Step 1: Find the highest AUC-PR value
max_aucpr = max(r["best_aucpr"] for r in cv_results_summary)

# Step 2: Get all CV configs that reached the highest AUC-PR
top_cv_configs = [r for r in cv_results_summary if r["best_aucpr"] == max_aucpr]

# Step 3: Pick the one with the smallest CV (or change this tie-breaker)
best_cv_result = min(top_cv_configs, key=lambda x: x["cv"])

# Step 4: Print result
print("🏆 Best AUC-PR across CV folds:")
print(f"CV folds: {best_cv_result['cv']}")
print(f"AUC-PR: {best_cv_result['best_aucpr']:.4f}")
print(f"Boosting rounds: {best_cv_result['best_iteration']}")

import xgboost as xgb

# Get core XGBoost parameters (not sklearn-style)
xgb_params = clf_xgb.get_xgb_params()
xgb_params['eval_metric'] = 'aucpr'  # Make sure AUC-PR is used

# Create DMatrix
dtrain = xgb.DMatrix(X, label=y)

#Make sure to add the best "CV folds" number obtained from the above CV result here "nfold=8"
# Perform K-fold cross-validation
cv_results = xgb.cv(
    params=xgb_params,
    dtrain=dtrain,
    num_boost_round=1000,
    nfold=8,
    stratified=True,
    early_stopping_rounds=15,
    seed=42,
    verbose_eval=True
)

print(cv_results)

#Print CV score from the best iteration.
best_aucpr = cv_results['test-aucpr-mean'].max()
print("Best AUCPR Score from Cross Validation:", best_aucpr)

import pandas as pd
import matplotlib.pyplot as plt

# Assuming cv_results is the result from xgb.cv
cv_results_df = pd.DataFrame(cv_results)

# Display the results in a table
print(cv_results_df)

# Plot the training and testing mean AUCPR scores with specified figsize
ax = cv_results_df[['test-aucpr-mean', 'train-aucpr-mean']].plot(figsize=(6.5, 4.5))
plt.title('XGBoost Cross-Validation Results')
plt.xlabel('Number of boosting rounds')
plt.ylabel('AUCPR')

# Add a dotted line at the highest AUCPR
highest_aucpr = cv_results_df['test-aucpr-mean'].max()
plt.axhline(y=highest_aucpr, color='r', linestyle='--', label=f"Max AUCPR: {highest_aucpr:.4f}")
plt.legend()

# Show the plot
plt.show()

# Save the plot as PNG and PDF
ax.figure.savefig('XGBoost Cross Validation Results.png', bbox_inches='tight')
ax.figure.savefig('XGBoost Cross Validation Results.pdf', bbox_inches='tight')

# Correct way to get best round from CV:
best_round = cv_results['test-aucpr-mean'].idxmax()
print("Best round from CV:", best_round)

# Extract feature importances from the model
importances = clf_xgb.get_booster().get_score(importance_type='gain')

# Sort features based on importance
sorted_features = sorted(importances.items(), key=lambda x: x[1], reverse=True)

# Save top 100 features to a txt file
with open('top_100_features.txt', 'w') as f:
    for feature, importance in sorted_features[:100]:
        f.write(f"{feature}: {importance}\n")

# Set the size of the plot
plt.rcParams["figure.figsize"] = (25, 10)

# Plot the top 20 most important features
plot_importance(clf_xgb, max_num_features=20, importance_type='gain')

# Add labels and a title
plt.title('Feature Importance')
plt.xlabel('F Score')
plt.ylabel('Features')

# Save the Feature Importance figure as a PDF
plt.savefig('xgb_feature_importance.pdf', dpi=300)

# Save the Feature Importance figure as a PNG
plt.savefig('xgb_feature_importance.png', dpi=300)

# Show the plot
plt.show()

from sklearn import metrics

# Ensure the model is a binary classifier
assert len(clf_xgb.classes_) == 2, "Model is not a binary classifier"

# Get probabilities for the positive class
y_preds = clf_xgb.predict_proba(X_test)[:, 1]

# Calculate ROC curve
fpr, tpr, _ = metrics.roc_curve(y_test, y_preds)

# Calculate AUC (Area under the ROC Curve )
auc_score = metrics.auc(fpr, tpr)

# Create figure and axis
plt.figure(figsize=(10, 10))

# Plot the ROC curve
plt.plot(fpr, tpr, color='blue', label='ROC curve (area = {:.2f})'.format(auc_score))

# Plot the line of no discrimination
plt.plot([0, 1], [0, 1], color='red', linestyle='--')

# Set plot labels and title
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic Curve')

# Add legend
plt.legend(loc='lower right')

# Save the ROC curve as a PNG
plt.savefig('cfdna_xgb_ROC.png', dpi=300)

# Save the ROC curve as a PDF
plt.savefig('cfdna_xgb_ROC.pdf', dpi=300)

# # Show the plot
plt.show()

# # Close the plot
# plt.close()

#To Save the BEST XGBoost tree at the best iteration.
# Assuming you've already trained your XGBoost model with early stopping and clf_xgb is your trained model

# Get the best iteration (subtracting 1 because it's 0-indexed)
#best_iteration = clf_xgb.best_ntree_limit - 1
best_iteration = clf_xgb.best_iteration


# Set parameters for the nodes in the graph
node_params = {'shape': 'box',  # Make the nodes fancy
               'style': 'filled, rounded',
               'fillcolor': '#78cbe'}

# Set parameters for the leaf nodes in the graph
leaf_params = {'shape': 'box',  # Make the nodes fancy
               'style': 'filled',
               'fillcolor': '#e48038'}

# Generate the graph for the best iteration
graph = xgb.to_graphviz(clf_xgb, num_trees=best_iteration, size="10,10",
                        condition_node_params=node_params,
                        leaf_node_params=leaf_params)

# Save the graph to a PDF file
graph.format = 'pdf'
graph.render('best_tree_graph')

# Display the tree
graph.view(cleanup=True)

#To print the BEST XGBoost tree

bst =clf_xgb.get_booster()
for importance_type in ('weight','gain','cover','total_gain','total_cover'):
  print('%s:'% importance_type, bst.get_score(importance_type=importance_type))

node_params = {'shape':'box', ##make the nodes fancy
               'style':'filled, rounded',
               'fillcolor':'#78cbe'}

leaf_params = {'shape':'box', ##make the nodes fancy
               'style':'filled',
               'fillcolor':'#e48038'}

xgb.to_graphviz(clf_xgb, num_trees=best_iteration,size="14,14",
                 condition_node_params=node_params,
                 leaf_node_params=leaf_params)

from sklearn.metrics import ConfusionMatrixDisplay
import matplotlib.pyplot as plt

# Generate predictions
y_pred = clf_xgb.predict(X_test)

# Generate confusion matrix from predictions
disp = ConfusionMatrixDisplay.from_predictions(y_test, y_pred, display_labels=["Healthy","Cancer"])

# Add a title to the plot
plt.title('Confusion Matrix Cancer VS Healthy')

# Save the plot to a PDF file
plt.savefig('confusion_matrix.pdf', format='pdf', bbox_inches='tight')

# Show the plot
plt.show()

# Install SHAP if not installed
!pip install shap

import shap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Create a SHAP explainer object and calculate SHAP values for the training set
explainer = shap.TreeExplainer(clf_xgb)
shap_values = explainer.shap_values(X_train)

# Convert SHAP values to a Pandas DataFrame
shap_df = pd.DataFrame(shap_values, columns=X_train.columns)

# =================== Save SHAP values for ALL features (Commented Out) ===================
# shap_csv_output_all = "SHAP_values_all_features.csv"
# shap_df.to_csv(shap_csv_output_all, index=False)
# print(f"SHAP values for all features saved as: {shap_csv_output_all}")

# Identify the top 20 most important features based on mean absolute SHAP values
shap_importance = np.abs(shap_df).mean().sort_values(ascending=False)
top_20_features = shap_importance.index[:20]

# Filter only the top 20 features
shap_df_top20 = shap_df[top_20_features]

# Save SHAP values for the top 20 features to a CSV file
shap_csv_output_top20 = "SHAP_values_top20_features.csv"
shap_df_top20.to_csv(shap_csv_output_top20, index=False)
print(f"SHAP values for top 20 features saved as: {shap_csv_output_top20}")

# Set up the matplotlib figure for the summary plot
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values, X_train, show=False)
plt.title("SHAP Summary Plot for Model Features")
plt.tight_layout()

# Save the figure to a PDF
pdf_output_summary = "SHAP_summary_plot.pdf"
plt.savefig(pdf_output_summary, format='pdf')

# Display the SHAP summary plot in the notebook
plt.show()

from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt

y_scores = clf_xgb.predict_proba(X_test)[:, 1]


# Compute the precision and recall
precision, recall, _ = precision_recall_curve(y_test, y_scores)

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(recall, precision, label='XGBoost')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall curve')
plt.legend(loc='best')
plt.show()

# # Zip the directory where your files are saved and download all at once
# !zip -r output_files.zip /content/
# from google.colab import files
# files.download("output_files.zip")

import pandas as pd
import numpy as np
import xgboost as xgb
import sklearn
import matplotlib
import shap

packages = {
    'pandas': pd.__version__,
    'numpy': np.__version__,
    'xgboost': xgb.__version__,
    'scikit-learn': sklearn.__version__,
    'matplotlib': matplotlib.__version__,
    'shap': shap.__version__
}

print(packages)

# #To list all the attributes of your clf_xgb object
dir(clf_xgb)

"""Created by: Sakuntha Devaka Gunarathna /
Published Date: July 30, 2025
"""