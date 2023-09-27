import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import math

def gaussian (x, FWHM, i , f):
    sigma = FWHM/2.355
    return np.exp(-np.power((x - i), 2.)/(2*sigma))*f

#FITTING A GAUSSIAN TO EACH LINE 
x = np.linspace(0.5, 6.0, 300)

archivo="SUMMARY_STEPS.dat"
out="spectrum_"+archivo
out2="spectrum_nm_"+archivo
data= np.loadtxt(archivo)
ro=data[:,2]
fo=data[:,3]

gaus = np.zeros(len(x))
for i in range(0,len(ro)):
    gaus += gaussian(x, 0.3, ro[i], fo[i])
with open(out, 'w') as f:
    f.write('#E (eV)    Intensity')
    f.write('\n')
    for i in range(len(x)):
        f.write(str(x[i]))
        f.write('   ')
        f.write(str(gaus[i]))
        f.write('\n')
with open(out2, 'w') as f:
    f.write('#E (nm)    Intensity')
    f.write('\n')
    for i in range(len(x)):
        f.write(str((10000000*27.211/219474)/x[i]))
        f.write('   ')
        f.write(str(gaus[i]))
        f.write('\n')
                                                                                                                                                                                    f.write('   ')
                                                                                                                                                                                                    f.write(str(gaus[i]))
                                                                                                                                                                                                                    f.write('\n')
                                                                                                                                                                                                                    
