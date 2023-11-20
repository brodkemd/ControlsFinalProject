import scipy.io

# Specify the MATLAB file path
mat_file_path = 'matrix.mat'

# Load the MATLAB file into a Python dictionary
mat_data = scipy.io.loadmat(mat_file_path)

# Access the variables stored in the MAT file
# For example, if you have a variable named 'matrix' in the MAT file
matrix = mat_data['k']

# Now you can use 'matrix' as a NumPy array in Python
print(matrix)