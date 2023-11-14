# FastICA in Python using sklearn
import numpy as np
from sklearn.decomposition import FastICA
import sys

# Load your data here, for example, a CSV file, or pass it directly.
# For demonstration purposes, I'm assuming the input is a CSV file.
data = np.loadtxt(sys.argv[1], delimiter=',')

# Initialize FastICA with number of components
ica = FastICA(n_components=3, random_state=0)

# Fit the model and transform the data
S_ = ica.fit_transform(data)  # Reconstructed signals
A_ = ica.mixing_  # Mixing matrix

# Save the independent components and the mixing matrix to CSV files
np.savetxt("independent_components.csv", S_, delimiter=',')
np.savetxt("mixing_matrix.csv", A_, delimiter=',')

# Optionally, print the mixing matrix to stdout to capture it in MATLAB
print(A_)
