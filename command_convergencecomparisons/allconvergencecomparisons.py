import math
import numpy as np
import subprocess

# 5 tests
#thetas_o = ["1.5"]
#alphas = ["0.6"]
#widths = ["0.0001", "0.0003", "0.0012", "0.005", "0.01"]

# In the paper:
# 45 tests
thetas_o = ["0.", "1.", "1.5"]
alphas = ["0.1", "0.25", "0.6"]
widths = ["0.0001","0.0003","0.0012","0.005","0.01"]

# n_runs = 8

# In the paper:
n_runs = 1000

pgf = 0

seed = 0

for width in widths:
    for alpha in alphas:
        for theta_o in thetas_o:
            print('Set of parameters ' + str(seed + 1))
            process = subprocess.run(['python', 'convergencecomparisons.py', theta_o, alpha, width, str(seed), str(
                n_runs), str(pgf)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
            print(process.stdout)
            seed = seed + 1
