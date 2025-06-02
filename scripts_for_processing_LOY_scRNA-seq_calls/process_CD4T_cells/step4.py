import pandas as pd
import phate
import matplotlib.pyplot as plt

# Step 1: Load Data
# Ensure that "data_matrix.csv" is in the same directory or provide the correct path.
print("Loading data matrix...")
data_matrix = pd.read_csv("data_matrix_CD4T_cells_for_phate.csv", index_col=0)  # Row names are treated as index

# Step 2: Perform PHATE
print("Performing PHATE analysis...")
#phate_op = phate.PHATE(n_components=2, knn=5, decay=40)  # Adjust parameters as needed
phate_operator = phate.PHATE(n_jobs=-2, n_pca=10)
phate_result = phate_operator.fit_transform(data_matrix)

# Step 3: Save PHATE Results
print("Saving PHATE results...")
phate_df = pd.DataFrame(phate_result, columns=["PHATE1", "PHATE2"], index=data_matrix.index)
phate_df.to_csv("phate_result_CD4T.csv", index=True)  # Save with row indices (e.g., cell names)

# Step 4: Plot PHATE Results
print("Plotting PHATE results...")
plt.figure(figsize=(8, 6))
plt.scatter(phate_df["PHATE1"], phate_df["PHATE2"], s=2, alpha=0.7, c='blue')
plt.title("PHATE Analysis")
plt.xlabel("PHATE1")
plt.ylabel("PHATE2")
plt.savefig("phate_plot_CD4T.png", dpi=300)  # Save plot as PNG
plt.show()

print("PHATE analysis and visualization completed.")

