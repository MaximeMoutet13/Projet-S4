__author__ = 'julesmichaud'
__filename__ = 'tests_jules.py'
__date__ = '26/05/20'

import numpy as np
import matplotlib.pyplot as plt

A = np.arange(1,10)
B = A[::-1]
C = np.concatenate(A,B)
print(A, B, C)
