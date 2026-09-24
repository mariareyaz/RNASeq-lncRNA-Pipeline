# =============================================================================
# Script 3: Machine Learning Framework for lncRNA Sequence Classification
# Description: Random Forest & Feature Importance modeling for distinguishing 
#              long non-coding RNAs from coding transcripts based on sequence features.
# =============================================================================

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score

# 1. Load Feature Matrix (ORF length, GC content, Fickett score, k-mer frequencies)
# Data extracted from candidate transcript FASTA sequences
df = pd.read_csv("transcript_features.csv")

# Features: ORF_len, GC_content, Fickett_score, 3mer_freq_1, 3mer_freq_2...
X = df.drop(columns=["transcript_id", "label"]) # label: 0 for lncRNA, 1 for mRNA
y = df["label"]

# 2. Train/Test Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# 3. Model Training (Random Forest Classifier)
clf = RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42)
clf.fit(X_train, y_train)

# 4. Evaluation
y_pred = clf.predict(X_test)
y_proba = clf.predict_proba(X_test)[:, 1]

print("--- AI Model Performance ---")
print(classification_report(y_test, y_pred))
print(f"ROC-AUC Score: {roc_auc_score(y_test, y_proba):.4f}")

# 5. Extract Feature Importance for Biomarker Discovery
feature_importances = pd.Series(clf.feature_importances_, index=X.columns).sort_values(ascending=False)
print("\nTop Predictors of Non-Coding Transcripts:")
print(feature_importances.head(5))
