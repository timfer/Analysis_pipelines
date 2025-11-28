#%%
import os
import numpy as np
import torch
import torch.nn as nn
from sklearn.model_selection import train_test_split

# 1. Load patient dictionaries


# 2. Inspect cohort and domain (batch) distributions
cohort_labels = np.array([d['cohort'] for d in patient_dicts])
batch_labels  = np.array([d['batch_id'] for d in patient_dicts])

unique_cohorts, cohort_counts = np.unique(cohort_labels, return_counts=True)
print("Cohorts:", dict(zip(unique_cohorts, cohort_counts)))

unique_batches, batch_counts = np.unique(batch_labels, return_counts=True)
print("Batches:", dict(zip(unique_batches, batch_counts)))

# 3. Split into train / validation at patient level
# Stratify on cohort to preserve class ratio
train_dicts, val_dicts = train_test_split(
    patient_dicts,
    test_size=0.2,
    stratify=cohort_labels,
    random_state=42
)

# 4. Compute class weights for cohort imbalance
total_patients = len(patient_dicts)
# counts per cohort from train set for more accurate weighting
train_cohorts = np.array([d['cohort'] for d in train_dicts])
unique_cohorts, train_counts = np.unique(train_cohorts, return_counts=True)
class_weights = total_patients / train_counts.astype(float)
# reorder to [weight_for_cohort0, weight_for_cohort1]
weight_tensor = torch.tensor(
    [class_weights[unique_cohorts.tolist().index(0)],
     class_weights[unique_cohorts.tolist().index(1)]],
    dtype=torch.float32
).to(device)

# 5. Initialize loss function with class weights
criterion = nn.CrossEntropyLoss(weight=weight_tensor)

# 6. Choose and adapt pretrained ResNet-18
from torchvision import models
resnet18 = models.resnet18(pretrained=True)
# Remove final FC and replace with identity (we'll add MIL head later)
resnet18.fc = nn.Identity()

# Wrap in your AttentionMIL model
model = AttentionMIL(backbone=resnet18, feature_dim=512, hidden_dim=128, num_cohorts=2)
model = model.to(device)

# 7. Summary of setup
print(f"Training on {len(train_dicts)} patients, validating on {len(val_dicts)} patients")
print("Class weights:", weight_tensor.cpu().numpy())
print("Model backbone: ResNet-18 (pretrained)")

# Next steps:
#  - Build DataLoader for train_dicts and val_dicts (patient-bag with dynamic HDF5 loading)
#  - Implement optimizer and training loop calling `criterion`
