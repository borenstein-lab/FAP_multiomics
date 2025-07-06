import numpy as np
import pandas as pd
import seaborn as sns
from plotnine import *
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from catboost import CatBoostClassifier, Pool
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import precision_recall_curve, average_precision_score, balanced_accuracy_score
from scipy.special import expit

def evaluate_model(model,model_name, X_test, y_test):
    """Evaluate model by calculating precision-recall values, PR AUC, and balanced accuracy."""
    if model_name == "SVM":
        y_probs = expit(model.decision_function(X_test))
    else:
        y_probs = model.predict_proba(X_test)[:, 1]
    precision, recall, _ = precision_recall_curve(y_test, y_probs)
    pr_auc = average_precision_score(y_test, y_probs)
    balanced_acc = balanced_accuracy_score(y_test, model.predict(X_test))
    return precision, recall, pr_auc, balanced_acc

def save_model_evaluation(results, model_name, precision, recall, pr_auc, balanced_acc):
    results[model_name]["precision"].append(precision)
    results[model_name]["recall"].append(recall)
    results[model_name]["pr_auc"].append(pr_auc)
    results[model_name]["balanced_acc"].append(balanced_acc)
    return results


def collect_feature_importances_df(X_columns, feature_importances):
    """Collect feature importances from all models and return as a reshaped DataFrame."""
    feature_importances_rows = []
    # Loop over models
    for model_name in feature_importances.keys():
        # Loop over each feature
        for i, feature in enumerate(X_columns):
            imp_mean = np.mean([imp[i] for imp in feature_importances[model_name]])
            imp_std = np.std([imp[i] for imp in feature_importances[model_name]])
            row = {"feature": feature,"model": model_name,"importance_mean": imp_mean,"importance_std": imp_std}
            feature_importances_rows.append(row)

    return pd.DataFrame(feature_importances_rows)


def evaluate_models(X, y, k=5, random_state=42):
    """Perform k-fold cross-validation for Logistic Regression, CatBoost, and SVM."""
    models = {
        "Logistic Regression": LogisticRegression(penalty="l1", solver="liblinear"),
        "SVM": LinearSVC(C=0.1, penalty="l1", dual=True),
        "CatBoost": None  # CatBoost is defined within the loop
    }

    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=random_state)

    # Store results
    results = {name: {"precision": [], "recall": [], "pr_auc": [], "balanced_acc": []} for name in models}
    feature_importances = {name: [] for name in models}
    
    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y), start=1):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # Compute class weights for the training fold
        class_weights = compute_class_weight(class_weight="balanced", classes=np.unique(y_train), y=y_train)
        class_weight_dict = {cls: weight for cls, weight in zip(np.unique(y_train), class_weights)}

        for model_name, model in models.items():
            # Train model (Catboost is a specieal case because we need to set class weights when defining the model)
            if model_name == "CatBoost":
                model = CatBoostClassifier(iterations=100,
                                           depth=3,
                                           learning_rate=0.1,
                                           loss_function="Logloss",
                                           scale_pos_weight=class_weights[1],
                                           l2_leaf_reg=10,
                                           verbose=False)
                train_pool = Pool(X_train, label=y_train)
                model.fit(train_pool)
            else:
                model.set_params(class_weight=class_weight_dict)
                model.fit(X_train, y_train)

            # Evaluate model and collect metrics
            precision, recall, pr_auc, balanced_acc = evaluate_model(model, model_name, X_test, y_test)
            results = save_model_evaluation(results, model_name, precision, recall, pr_auc, balanced_acc)

            # Extract and store feature importances
            if hasattr(model, 'coef_'):  # Logistic regression (coefficients)
                feature_importances[model_name].append(np.abs(model.coef_[0]))
            elif hasattr(model, 'feature_importances_'):  # CatBoost (feature importances)
                feature_importances[model_name].append(model.feature_importances_)

    # Summarize feature importances into a DataFrame
    feature_importances_df = collect_feature_importances_df(X.columns, feature_importances)

    # Plot PR curves for all models
    #plot_pr_curves(results)

    return results, feature_importances_df

def evaluate_models_n_times(X, y, n_splits, n_times):
    """
    Evaluate models multiple times with k-fold cross-validation.
    
    Parameters:
    X (DataFrame): Feature data.
    y (Series): Target variable.
    n_splits (int): Number of splits for cross-validation.
    n_times (int): Number of times to repeat the evaluation.

    Returns:
    DataFrame: Model performance statistics.
    DataFrame: Feature importance statistics for stable features within a seed.
    """
    # Initialize empty DataFrames to store results
    pr_df_all = pd.DataFrame()
    feat_imp_all = pd.DataFrame()
    for seed in [42 + x for x in range(n_times)]:
        # Get model performance statistics and feature importance for k-fold cross-validation for the current seed
        results, feat_imp = evaluate_models(X,y, n_splits,seed)
        # Filter feature importance to only include features with mean importance greater than std
        feat_imp = feat_imp[feat_imp.importance_mean > feat_imp.importance_std]
        feat_imp['seed'] = f"Seed {seed}"
        feat_imp_all = pd.concat([feat_imp_all, feat_imp], ignore_index=True)
        
        # Calculate mean and std statistics across the folds of each model for the current seed
        for model, model_res in results.items():
            precision_mean = np.mean(np.array(model_res['precision']), axis=0)
            recall_mean = np.mean(np.array(model_res['recall']), axis=0)
            pr_auc = np.mean(model_res['pr_auc'])
            pr_auc_std = np.std(model_res['pr_auc'])
            balanced_acc = np.mean(model_res['balanced_acc'])
            balanced_acc_std = np.std(model_res['balanced_acc'])
            
            pr_df = pd.DataFrame({'ind': range(len(precision_mean)),'model': model,
                                'precision': precision_mean,'recall': recall_mean,
                                 'pr_auc': pr_auc, 'pr_auc_std': pr_auc_std,
                                 'balanced_acc': balanced_acc, 'balanced_acc_std': balanced_acc_std})
            pr_df['seed'] = f"Seed {seed} (PR AUC = {pr_auc:.2f}±{pr_auc_std:.2f})"
            pr_df_all = pd.concat([pr_df_all, pr_df], ignore_index=True)
            
    return(pr_df_all, feat_imp_all)

def plot_avg_pr_curves(df): 
    models = df["model"].unique()
    n_models = len(models)
    fig, axes = plt.subplots(1, n_models, figsize=(6 * n_models, 5), sharey=True)

    if n_models == 1:
        axes = [axes]  # Ensure axes is iterable even if one subplot

    # Loop through each model to create its subplot
    for ax, model in zip(axes, models):
        df_model = df[df["model"] == model]
        aucs = df_model['pr_auc'].drop_duplicates()

        # Plot one line per seed with confidence interval
        sns.lineplot(data=df_model, x="recall",y="precision",
            hue="seed", errorbar=None, linewidth=4, alpha=0.5,palette=sns.color_palette("Set1", max(n_models,5)),  ax=ax)

        ax.set_title(f"{model}: Mean PR AUC = {np.mean(aucs):.2f}±{np.std(aucs):.2f}")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.legend(title="Seed", loc="lower left", bbox_to_anchor=(0, 0))
        

    plt.tight_layout()
    plt.show()

## to be deleted
def plot_pr_curves(results):
    """Plot Precision-Recall curves for each model."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (model_name, data) in zip(axes, results.items()):
        for fold_idx, (precision, recall, pr_auc) in enumerate(zip(data["precision"], data["recall"], data["pr_auc"]), start=1):
            ax.plot(recall, precision, label=f"Fold {fold_idx} (PR AUC={pr_auc:.2f})")

        ax.set_title(f"{model_name} - PR Curve")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.legend()
        ax.grid(True)

    plt.tight_layout()
    plt.show()