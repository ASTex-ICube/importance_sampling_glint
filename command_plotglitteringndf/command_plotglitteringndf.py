import math
import numpy as np
import subprocess

imagesize = 256

alphax = 0.6
alphay = 0.6

# 147 microfacets
width = 0.00052
# 2,040 microfacets
# width = 0.002
# 161,686 microfacets
# width = 0.017

outfilename = "plotglitteringndf.py"

process = subprocess.run(['glinttool', 'plotglitteringndf', '-imagesize', str(imagesize), '-alphax', str(alphax), '-alphay', str(alphay), '-dsdx', str(width),'-dtdy', str(width), 'dict_N768_nLevels8.exr',
                          outfilename], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
print(process.stdout)
process = subprocess.run(['python', outfilename],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
print(process.stdout)
