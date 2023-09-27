import tkinter as tk
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from modules.gaussians import Gaussians
from modules.plot_gaussians_gui import Gui_app


def ev_to_nm(E):
    return 1240./E

def parse_spectrum_orca(infile, units = "nm"):
    """Read Orca output file. Return spectrum. If units is not nanometers,
       the units are transformed on the fly.
    """

    with open(infile, "r") as f:
        lines = f.read().splitlines()
    E, f = [], []
    for cn1, iline in enumerate(lines):
        if("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in iline):
            for cn2, jline in enumerate(lines[cn1 + 5:]):
                brk_cond = ("----------" in jline) or (len(jline) == 0)
                if(brk_cond):
                    break
#                print(jline)
                row = jline.split()
                E.append(float(row[2]))
                f.append(float(row[3]))
#    return(np.array(E), np.array(f))
    if(units.lower() == "ev"):
        E = [ev_to_nm(i) for i in E]
    return(E, f)


if(__name__ == "__main__"):
    from argparse import ArgumentParser

    # Set DPI for matplotlib
#    plt.rcParams["figure.dpi"] = 150
    parser = ArgumentParser(
           "Plot spectrum from an ensemble of Orca calculations.\n\
            Retrieve each output from the output geomN/geomN.inp.out,\n\
            where N is the geometry of interest (A range can also \n\
                    be defined)")

    parser.set_defaults(units = "eV", color_bars = False, max = 5)    
    parser.add_argument("-f", nargs="+", help = "Range of geometries to analyze")
    parser.add_argument("-u", dest = "units", type = str,\
            help = "Units of energy (default = eV)")
    parser.add_argument("-m", dest = "max", type = int,\
            help = "maximum state to consider (default = 5)")
    parser.add_argument("-cb", dest = "color_bars", type = bool,\
            help = "Color bars by state (default = False)")

    options = parser.parse_args()
    rng     = options.f   
    if(len(rng) == 3):
        start, stop, step = [int(i) for i in  rng]
        stop += 1
    elif(len(rng) == 2):
        start, stop = [int(i) for i in rng]
        step = 1
        stop += 1
    elif(len(rng) == 1):
        start = int(rng[0])
        stop  = start + 1
        step  = 1

    color_bars  = options.color_bars
    units       = options.units
    geoms       = list(range(start, stop, step))
    max_state   = options.max
    
    # Define x and y arrays
    # that is, energy and f.
    # By default, E is in nm. If units 
    # are eV, units are transformed on initialization
    x, y = [], []
    for cnt, igeom in enumerate(geoms):
        infile = "geom{:d}/geom{:d}.inp.out".format(igeom, igeom)
        xi, yi = parse_spectrum_orca(infile, units = units)
        x.append(xi[:max_state])
        y.append(yi[:max_state])
    
#    x = np.array(x, dtype = object)
#    y = np.array(y, dtype = object)
    x = np.array(x)
    y = np.array(y)
    print(x)
    print(type(x))


    root = tk.Tk()
    app = Gui_app(root, x, y, units = units, color_bars = color_bars, max_state = max_state) # Instantiate App
    root.mainloop()



