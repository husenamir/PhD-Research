# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 21:33:14 2023

@author: ahusen
"""

import numpy as np

def conjugate_gradient(S, b):
    N = b.shape[0]
    x = np.zeros(b.shape)
    r = b
    d = r
    for k in range(1, N+1):
        alpha = np.matmul(r.transpose(), r)/np.matmul(d.transpose(), np.matmul(S, d))
        x = x + alpha * d
        r_next = r - alpha * np.matmul(S, d)
        u = np.matmul(r_next.transpose(), r_next)
        d = np.matmul(r.transpose(), r)
        beta = u/d
        d = r_next + beta * d
        r = r_next
    return x

A = np.array([[2,-1,0],[-1,2,-1], [0,-1,2]])
b = np.array([[1],[0],[0]])

x = conjugate_gradient(A, b)
print('After ', len(b), 'iterations, Solution: ', x)
