import numpy as np

class Gaussians(object):
    """Convolution of a bar spectrum with gaussians. 
       Good values for xo and xf could be xo = min(x) - 3*FWHM
       and xf = max(x) + 3*FWHM"""
    def __init__(self, barsx, barsy, xo, xf, mesh, FWHM):
        self.sigma = FWHM/2.355
        self.barsx = barsx
        self.barsy = barsy
        self.x = np.arange(xo, xf, 1./mesh)
        self.y = np.zeros(self.x.size)

    def gen_gaussians(self):
        """Self explanatory"""
        for E, f in zip(self.barsx, self.barsy):
#            print(E, f)
            center_gauss = np.full(self.x.size, E)
            gauss = f*np.exp(-np.power(self.x - center_gauss, 2)/(2.*self.sigma**2))
            self.y += gauss
            gauss = np.zeros(self.x.size)


