import math
import numpy as np
import subprocess
import sys

###############################################################################
# Simple function extracting a string from a number to be included in a file name


def toStringForFilename(number):
    out_string = ''
    if number < 0:
        out_string += 'm'
    number_str = str(number)
    if number_str.find('e') != -1:
        out_string += str(number)
    else:
        out_string += str(int(number))
        decimal_part = abs(number) - math.floor(abs(number))
        if decimal_part == 0:
            return out_string
        decimal_part_str = str(decimal_part)
        out_string += 'p'
        out_string += decimal_part_str[2:6]
    return out_string
###############################################################################


###############################################################################
# Processes command arguments
if len(sys.argv) > 1:
    theta_o = float(sys.argv[1])
else:
    theta_o = 1.5

if len(sys.argv) > 2:
    alpha = float(sys.argv[2])
else:
    alpha = 0.6

if len(sys.argv) > 3:
    pixel_width = float(sys.argv[3])
else:
    pixel_width = 0.0001

if len(sys.argv) > 4:
    stx = float(sys.argv[4])
else:
    stx = 0.

if len(sys.argv) > 5:
    n_runs = int(sys.argv[5])
else:
    n_runs = 8

if len(sys.argv) > 6:
    pgf = int(sys.argv[6])
else:
    pgf = 0
###############################################################################

# File name without extension
outfilename_without_extension = 'theta_o' + toStringForFilename(theta_o)
outfilename_without_extension += '_alpha' + toStringForFilename(alpha)
outfilename_without_extension += '_pw' + toStringForFilename(pixel_width)

# File name with extension
outfilename = outfilename_without_extension + '.py'

# Call up the glitter tool
# glinttool must be in the path
# The dictionary dict_N768_nLevels8.exr must be in the python script folder
print('outfilename=' + outfilename)
process = subprocess.run(['glinttool', 'convergencecomparisons', '-stx', str(stx), '-thetao', str(theta_o), '-alphax', str(alpha), '-alphay', str(alpha), '-dstdxx', str(pixel_width), '-dstdyy', str(pixel_width), '-pgf', str(pgf), '-nruns', str(n_runs),
                          'dict_N768_nLevels8.exr', outfilename], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
print(process.stdout)

# Matplotlib: radiance as a function of the number of samples
process = subprocess.run(['python', outfilename_without_extension + '_radiance_n.py'],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
print(process.stdout)

# Matplotlib: pointwise boxplots
process = subprocess.run(['python', outfilename_without_extension + '_pointwise_boxplots.py'],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
print(process.stdout)