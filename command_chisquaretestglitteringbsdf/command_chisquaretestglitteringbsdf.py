import math
import numpy as np
import subprocess

# BSDFs need more samples and a bigger integration resolution
nsamplehisto = 256000000
res = 4096

# Not enough samples and too low resolution. Use these parameters only to see the BSDF plot
# nsamplehisto = 1000000
# res = 512

alphax = 0.25
alphay = 0.25
rho = 0.

stx = 0.

width = 0.001

thetao = 0.2
phio = 0.

# microfacet relative area
mra = 1.0

outfilenameWithoutExtension = "plotchisquaretestglitteringbsdf"
outfilename = outfilenameWithoutExtension + ".py"

process = subprocess.run(['glinttool', 'chisquaretestglitteringbsdf', '-nsamplehisto', str(nsamplehisto), '-res', str(res), '-alphax', str(alphax), '-alphay', str(alphay), '-rho', str(rho), '-stx', str(stx), '-dstdxx', str(width), '-dstdyy', str(width), '-thetao', str(thetao), '-phio', str(phio), '-mra', str(mra), 'dict_N768_nLevels8.exr',
                          outfilename], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
print(process.stdout)

# Memory error with res = 4096. Too heavy for matplotlib. 
# process = subprocess.run(['python', outfilenameWithoutExtension + "_test1.py"],
#                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
# print(process.stdout)
