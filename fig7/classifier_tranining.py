import pandas as pd
from mrmr import mrmr_classif
import xgboost as xgb
import numpy as np
from sklearn.model_selection import train_test_split,cross_val_score,StratifiedKFold
from sklearn.metrics import accuracy_score,roc_auc_score
from sklearn.utils.class_weight import compute_sample_weight
import joblib
def feture_selection():
    list_genes = []
    fr = open("gene_overlap.tsv","r")
    for line in fr.readlines():
        gene = line.strip()
        list_genes.append(gene)

    dict_cluster = {}
    list_id = []
    fr0 = open("clustering_results.tsv","r")
    for line in fr0.readlines():
        if line.strip().startswith("ID"):
            continue
        else:
            sp = line.strip().split("\t")
            list_id.append(sp[1])
            dict_cluster[sp[1]] = int(sp[2])-1

    df = pd.read_csv("CPTAC_cpm_filtered.tsv",delimiter="\t",index_col=0)
    df = df[list_id]
    df_transposed = df.T
    print(df_transposed)
    filtered_df = df_transposed[list_genes]

    y = []
    index = []
    for id in filtered_df.index:
        y.append(dict_cluster[id])
        index.append(id)

    y_series = pd.Series(y,index=index)

    list_k = [10,20,50,100,200,300,400,500,1000,2000,3000]
    fw0 = open("performance_diff_k.tsv","w")
    fw0.write("Num\tAUC\n")
    rng = np.random.default_rng()
    for iteration in range(10):
        print(iteration)
        random_seed = rng.integers(low=0, high=1e9)
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_seed)
        for train_index, test_index in skf.split(filtered_df, y_series):
            X_train, X_test = filtered_df.iloc[train_index], filtered_df.iloc[test_index]
            Y_train, Y_test = y_series.iloc[train_index], y_series.iloc[test_index]
            X_train = X_train.apply(zscore)
            X_test = X_test.apply(zscore)
            top_features = mrmr_classif(X_train, Y_train, 3000)
            print(top_features[0:10])
            for i in range(len(list_k)):
                k = list_k[i]
                # Feature selection here with k features
                selected_features = top_features[:k]

                # Select the features from train and test set
                X_train_fs = X_train[selected_features]
                X_test_fs = X_test[selected_features]

                # Compute sample weights
                sample_weights = compute_sample_weight(class_weight="balanced", y=Y_train)

                # Initialize and train the XGBoost classifier
                model = xgb.XGBClassifier(objective='multi:softprob', num_class=8)
                model.fit(X_train_fs, Y_train, sample_weight=sample_weights)

                # Predict and evaluate the model
                Y_pred = model.predict_proba(X_test_fs)
                auc = roc_auc_score(Y_test, Y_pred, average="weighted", multi_class="ovr")
                print(str(k) + "\t" + str(auc))
                fw0.write(str(k) + "\t" + str(auc) + "\n")

    fw0.flush()
    fw0.close()

def model_traing(k):
    list_genes = []
    fr = open("gene_overlap.tsv", "r")
    for line in fr.readlines():
        gene = line.strip()
        list_genes.append(gene)

    dict_cluster = {}
    list_id = []
    fr0 = open("clustering_results.tsv", "r")
    for line in fr0.readlines():
        if line.strip().startswith("ID"):
            continue
        else:
            sp = line.strip().split("\t")
            list_id.append(sp[1])
            dict_cluster[sp[1]] = int(sp[2]) - 1

    df = pd.read_csv("Michigan_cpm_filtered.tsv",delimiter="\t", index_col=0)
    df = df[list_id]
    df_transposed = df.T
    print(df_transposed)
    filtered_df = df_transposed[list_genes]
    filtered_df = filtered_df.apply(zscore)
    print(filtered_df)
    y = []
    index = []
    for id in filtered_df.index:
        y.append(dict_cluster[id])
        index.append(id)

    y_series = pd.Series(y, index=index)

    selected_features = mrmr_classif(filtered_df,y_series,k)
    print(selected_features)

    df_training = filtered_df[selected_features]
    print(df_training)
    fw = open("feature_mrmr_" + str(k) + ".tsv", "w")
    for i in range(len(df_training.columns)):
        fw.write(df_training.columns[i] + "\t" + str(i+1) + "\n")
    fw.flush()
    fw.close()

    sample_weights = compute_sample_weight(class_weight="balanced",y=y)
    # Initialize the XGBoost classifier
    model = xgb.XGBClassifier(objective='multi:softprob',num_class=8)
    # Train the model with the entire training set
    model.fit(df_training, y,sample_weight=sample_weights)

    # Save the model
    joblib.dump(model, "best_model_" + str(k) + ".pkl")

def preidct(k):
    dict_genes = {}
    selected_genes = []
    fr = open("feature_mrmr_" + str(k) + ".tsv", "r")
    for line in fr.readlines():
        sp = line.strip().split("\t")
        gene = sp[0]
        order = int(sp[1])
        dict_genes[gene] = order
        selected_genes.append(gene)

    df = pd.read_csv("BeatAML_cpm_filtered.tsv",delimiter="\t", index_col=0)
    df_transposed = df.T
    filtered_df = df_transposed[selected_genes]
    sorted_columns = sorted(filtered_df.columns, key=lambda x: dict_genes[x])
    filtered_df = filtered_df[sorted_columns]
    filtered_df = filtered_df.apply(zscore)
    list_id = []
    for id in filtered_df.index:
        list_id.append(id)
    model = joblib.load("best_model_" + str(k) + ".pkl")

    Y_pred = model.predict(filtered_df)

    fw = open("Prediction_results_" + str(k) + ".tsv","w")
    fw.write("ID\tPred_cluster\n")
    for i in range(len(list_id)):
        fw.write(list_id[i] + "\t" + str(Y_pred[i]+1) + "\n")
    fw.flush()
    fw.close()

