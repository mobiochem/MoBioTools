import numpy as np
import sys
infile = sys.argv[1]
data = np.genfromtxt(infile, usecols = (0), skip_header = 1, unpack = True)
num = len(np.unique(data))
print(num)
