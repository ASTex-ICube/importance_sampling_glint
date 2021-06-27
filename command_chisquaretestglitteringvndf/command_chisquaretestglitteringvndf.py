import math
import numpy as np
import subprocess

nsamplehisto = 1000000
res = 512

alphax = 0.3
alphay = 0.3
rho = 0.

stx = 0.

width = 0.001

thetao = 1.5
phio = 0.

outfilenameWithoutExtension = "plotchisquaretestglitteringvndf"
outfilename = outfilenameWithoutExtension + ".py"

process = subprocess.run(['glinttool', 'chisquaretestglitteringvndf', '-nsamplehisto', str(nsamplehisto), '-res', str(res), '-alphax', str(alphax), '-alphay', str(alphay), '-rho', str(rho), '-stx', str(stx), '-dstdxx', str(width), '-dstdyy', str(width), '-thetao', str(thetao), '-phio', str(phio), 'dict_N768_nLevels8.exr',
                          outfilename], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
print(process.stdout)
process = subprocess.run(['python', outfilenameWithoutExtension + "_test1.py"],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
print(process.stdout)
