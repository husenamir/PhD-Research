#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:10:55 2025

@author: amirhusen
"""

import numpy as np
from matplotlib import pyplot as plt

"""
excellent_a = np.array([3.0, 2.8,  2.6, 2.8, 3.2, 3.0, 2.6, 2.6, 2.8, 3.0])
excellent_m = np.array([50, 50, 60,  70,  70,  80,  90,  100, 100, 100])

good_a = np.array([2.6, 3.2, 2.8, 3.2, 2.6, 2.6, 2.8, 3.2])
good_m = np.array([50,  50,  60,  60,  70,  80,  80,  100])

bad_a = np.array([3.0, 3.0, 3.2, 2.8, 3.0, 3.2])
bad_m = np.array([60,  70,  80,  90,  90,  90])

reference_a = np.array([3.03,  2.88,   2.86, 3.15])
reference_m = np.array([50.94, 51.99,  55.85, 95.94])


fig, ax = plt.subplots()
ax.plot(excellent_m, excellent_a,  label = 'Stable', linestyle = 'None', marker = 'P', color = 'green')
ax.plot(good_m, good_a,  label = 'Almost Stable', linestyle = 'None', marker = '^', color = 'blue')
ax.plot(bad_m, bad_a, label = 'Unstable', linestyle = 'None', marker = 'X', color = 'red')
ax.plot(reference_m, reference_a,  label = 'Reference', linestyle = 'None', marker = 'o', color = 'black')

ax.set_xlabel('Mass')
ax.set_ylabel('Lattice Parameter')
ax.legend()
"""
import re

file_path = "/Users/amirhusen/Desktop/Amir/Research/phonopy/2_6_50/generation_output.txt"

with open(file_path, "r") as f:
    content = f.read()

# Extract all "Best solution" blocks
matches = re.findall(r"Best solution:\s*\[([^\]]+)\]", content, re.DOTALL)
last_solution = matches[-1] if matches else ""

# Convert text → float values
numbers = [float(x) for x in last_solution.replace("\n", " ").split()]

# Format to 8 decimal places
# 1. Comma-separated string
comma_separated = ", ".join([f"{num:.8f}" for num in numbers])

# 2. Python list (floats, not strings)
python_list = [float(f"{num:.8f}") for num in numbers]

print(numbers)
print(comma_separated)
print("Python list:", python_list)
