#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 09:47:05 2023
First I will make self interacting random walk in 1D, 2D and 3D and try to 
reproduce the results of this guy. I
I will work in discrete lattice for simplicity. 
Then I will calculate the First Passage times and other things. 

How much sharp jumps? Enough that at time less than Tjump, behaviour is normal. 

Do I start with confined or infinite geometry? 
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the file into a pandas DataFrame
df = pd.read_csv('first_passage_times.txt', names=['First Passage Times'])

# Plot the distribution
plt.figure(figsize=(10, 6))
plt.hist(df['First Passage Times'], bins=np.logspace(0,5,100), color='skyblue', edgecolor='black', alpha=0.7)
plt.title('Distribution of First Passage Times')
plt.xlabel('First Passage Time')
plt.ylabel('Frequency')
plt.yscale('log')
plt.xscale('log')
plt.grid(True)
plt.show()
#%%
import numpy as np
import matplotlib.pyplot as plt

# Load data
first_passage_times = np.loadtxt("first_passage_times.txt")

# Total number of simulations
total_simulations = len(first_passage_times)

# Define the time range
max_time = np.max(first_passage_times)
times = np.logspace(0, np.log10(max_time), 1000)

# Calculate survival probabilities
survival_probabilities = []
for t in times:  
    survival_probabilities.append(np.sum(first_passage_times > t) / total_simulations)

# Plot
plt.figure(figsize=(8, 6))
plt.plot(times, survival_probabilities)
plt.xscale('log')
plt.xlabel('Time')
plt.ylabel('Survival Probability)
plt.grid(True)
plt.yscale('log')

plt.title('Survival Probability in 2D for $beta =1$ vs Time(t)')
plt.show()

#%%


# Replace with your actual file name
filename = "Sim10.txt"

# Read the data
with open(filename, 'r') as f:
    data = f.readlines()

# Split the data into x and y coordinates
x = []
y = []
for line in data:
    xi, yi = line.split()
    x.append(float(xi))
    y.append(float(yi))

# Plot the data
plt.plot(x, y)
plt.title('Trajectory Plot')
plt.xlabel('X')
plt.ylabel('Y')
# plt.xlim(0, 1000) 
# plt.ylim(0, 1000) 
plt.show()#%% Plot Sparse Matrix2
#%%
from scipy.sparse import coo_matrix

# Load the data from the text file
data = np.loadtxt('Gsim1.txt')

# Get the rows, cols, and values from the data
rows = data[:,0].astype(int)
cols = data[:,1].astype(int)
values = data[:,2]

# Define the shape of the matrix
shape = (100, 100)

# Create a sparse matrix from the rows, cols, and values
sparse_matrix = coo_matrix((values, (rows, cols)), shape=shape)

# Convert the sparse matrix to a dense matrix for plotting
dense_matrix = sparse_matrix.toarray()

# Create the plot
plt.imshow(dense_matrix, cmap='hot', interpolation='nearest')
plt.show()



