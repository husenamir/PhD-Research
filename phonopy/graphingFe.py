#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 16:59:17 2025

@author: amirhusen
"""


"""
data = pd.read_csv('/Users/amirhusen/Downloads/Default Dataset.csv', header=None, names= ['x', 'y'])
x_data1 = data['x']
y_data1 = data['y']

# Load the CSV data
data = pd.read_csv('/Users/amirhusen/Downloads/Default Dataset.csv', header=None)
x = data[0]
y = data[1]


# Plotting the data
plt.figure(figsize=(10, 5))
plt.plot(x, y, marker='o', linestyle='', color='b')
plt.title('Plot of y_data1 over x_data1')
plt.xlabel('Index (or some unit of X)')
plt.ylabel('Values of y_data1')

plt.show()
"""
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import fitz  # PyMuPDF
import pandas as pd

# Load the PDF and extract the first page as an image
pdf_path = "/Users/amirhusen/Desktop/band_ga.pdf"
doc = fitz.open(pdf_path)
page = doc[0]
pix = page.get_pixmap()
img = Image.open(io.BytesIO(pix.tobytes("png")))

# Load the CSV dataset
df = pd.read_csv('/Users/amirhusen/Downloads/Default Dataset.csv', header=None)
x_values = df.iloc[:, 0]  # First column as x
y_values = df.iloc[:, 1]  # Second column as y

# Determine the scale
x_min, x_max = min(x_values), max(x_values)
y_min, y_max = min(y_values), max(y_values)

# Plot the existing graph as a background
plt.figure(figsize=(8, 6))
plt.imshow(img, extent=[x_min, x_max, y_min, y_max], aspect='auto')

# Overlay new dataset
plt.scatter(x_values, y_values, color='blue', marker='o', label='New Data')

# Labels and legend
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.legend()
plt.title("Overlay of New Data on Existing Graph")

# Remove grid
plt.grid(False)

# Show plot
plt.show()
