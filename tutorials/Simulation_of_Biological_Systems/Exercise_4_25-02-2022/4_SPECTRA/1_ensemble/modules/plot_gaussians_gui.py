import tkinter as tk
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,\
        NavigationToolbar2Tk
from .gaussians import Gaussians

class Gui_app(object):
    
    def __init__(self, root, x, y, units = "ev", color_bars = False, max_state = 5):
        """Args:
           root = tk.Tk() object
           x    = array of energies with shape (N,S)
                  N = number of geometries
                  S = number of states

           y    = array of fosc with shape (N,S)
                  N = number of geometries
                  S = number of states
        """
        self.ngeoms  = len(x)
#        self.nstates = len(x[0])
        self.nstates = max_state
##        self.nstates = len(x)
        self.E = x[:,:max_state].flatten()
        self.f = y[:,:max_state].flatten()

        self.E_bar   = np.transpose(x)[:max_state]
        self.f_bar   = np.transpose(y)[:max_state]
        self.frame   = tk.Frame(root)
        self.wentry  = tk.Entry(self.frame) # fwhm entry
        self.wentry.bind("<Return>", self.get_fwhm_entry) # bind to gaussian plot command
        self.wentry.pack(side = "left")
        self.wbutton = tk.Button(self.frame, text = "Get FWHM",\
                                 command = self.get_fwhm)
        self.wbutton.pack(side = "right")
        
        # Other input arguments to vary the spectrum
        self.color_bars = color_bars
        self.units = units.lower()

        # Instantiate figure and canvas
        # General ax settings
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlabel("E [{:s}]".format(units))
        self.ax.set_ylabel("f")
        self.canvas = FigureCanvasTkAgg(self.fig, master = root)
        self.canvas.get_tk_widget().pack(side = "top", fill = "both", expand = 1)

        toolbar = NavigationToolbar2Tk(self.canvas, root)
        toolbar.update()
        self.canvas.draw()
        self.frame.pack()

        self.plot_bars()

    def get_fwhm(self):
        try:
            FWHM = float(self.wentry.get())
        except:
            raise ValueError("Insert a valid FWHM")
        xo   = np.min(self.E) - 3*FWHM 
        xf   = np.max(self.E) + 3*FWHM 
        gau = Gaussians(self.E, self.f, xo, xf, 200, FWHM)
        gau.gen_gaussians()
        # Delete previous gaussians
        try:
            self.ax.lines[0].remove()
        except:
            pass
        self.ax.plot(gau.x, gau.y)
        self.canvas.draw()

    def get_fwhm_entry(self, event):
        try:
            FWHM = float(self.wentry.get())
        except:
            raise ValueError("Insert a valid FWHM")
        xo   = np.min(self.E) - 3*FWHM 
        xf   = np.max(self.E) + 3*FWHM 
        gau = Gaussians(self.E, self.f, xo, xf, 200, FWHM)
        gau.gen_gaussians()
        # Delete previous gaussians
        try:
            self.ax.lines[0].remove()
        except:
            print("No gaussians to be deleted.")
            pass
        self.ax.plot(gau.x, gau.y, color = "gray")
        self.canvas.draw()


    def plot_bars(self):
        widths = {"ev": 0.01,
                  "nm": 2.0}
        w = widths[self.units]
        print("WIDTH = ", w)
        if(self.color_bars):
            colors = [cm.gnuplot(i) for i in np.linspace(0,1,self.nstates)]
            for i in range(0, self.nstates):
                x = self.E_bar[i]
                y = self.f_bar[i]
                color = colors[i]
                self.ax.bar(x, y, color = color, width = w, label = "S{}".format(i+1))
        else:
            self.ax.bar(self.E, self.f, color = "black", width = w)
            
#        self.ax.bar(self.E, self.f, color = "black", width = 0.01)
#        print(self.E)
#        self.ax.bar(self.E, self.f, color = "black", width = w)
        if(self.units == "ev"):
            self.ax.legend(loc = "upper left")
        elif(self.units == "nm"):
            self.ax.legend(loc = "upper right")
        else:
            print("No legend is being printed.")

        self.canvas.draw()
        



if(__name__ == "__main__"):
    from argparse import ArgumentParser
    parser = ArgumentParser("GUI for dynamically changing the broadness of the gaussians of a spectrum")
    parser.set_defaults(units = "eV", color_bars = False, max = 5)    
    parser.add_argument("-f", dest = "infile", type = str,\
            help = "Input file on xy format")
    parser.add_argument("-u", dest = "units", type = str,\
            help = "Units of energy (default = eV)")
    parser.add_argument("-m", dest = "max", type = int,\
            help = "maximum state to consider (default = 5)")
    parser.add_argument("-cb", dest = "color_bars", type = bool,\
            help = "Color bars by state (default = False)")

    options     = parser.parse_args()
    infile      = options.infile
    units       = options.units
    color_bars  = options.color_bars
    max_state   = options.max
    print(color_bars)

    x, y    = np.loadtxt(infile, usecols = (0,1), unpack = True)
    x = np.reshape(x, (1,len(x)))
    y = np.reshape(y, (1,len(y)))

    root = tk.Tk()
    app = Gui_app(root, x, y, units = units, color_bars = color_bars, max_state = max_state) # Instantiate App
    root.mainloop()


