import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

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
    import os
    # Plot "convergence" of each plot as a scatter plot
    # energy in eV
    currdir   = os.getcwd() + "/"
    max_state = 6
    igeom     = 6
    units     = "eV"
    
    # Consider up to 5 residues
    res     = list(range(1,6))
    states  = list(range(1,7))
    nstates = len(states)

    # Parse energies
    x = []

    for ires in res:
        infile = "{:d}res-geom{:d}/geom{:d}.inp.out".format(ires,igeom, igeom)
        xi, yi = parse_spectrum_orca(infile, units = units)
        x.append(xi[:max_state])

    x = np.transpose(x)

    # Initialize plot
    plt.ion()
    fig, ax = plt.subplots()
    ax.set_xlabel("Number of water molecules")
    ax.set_ylabel("E [eV]")
    colors = [cm.gnuplot(i) for i in np.linspace(0, 1, nstates)]

    for cnt, istate in enumerate(states):
        ax.scatter(res, x[cnt], color = colors[cnt], label = "S{:d}".format(istate))

    ax.legend(loc = "upper right")
    fig.savefig("comparison.png", dpi = 500)







    

