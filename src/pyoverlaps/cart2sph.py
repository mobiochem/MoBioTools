#!/usr/bin/env python3
from math import pi

#pi = 3.141592653589793;
offcrt = [0, 1, 10, 40]
trans_cart2sph = \
        [
        #S
        (1./(4*pi))**0.5,
        #PX
        (3./(4*pi))**0.5,
        0.,
        0.,
        #PY
        0.,
        (3./(4*pi))**0.5,
        0.,
        #PZ
        0.,
        0.,
        (3./(4*pi))**0.5,
        #DXY
        0.,
        1.092548430592079070,
        0.,
        0.,
        0.,
        0.,
        #DYZ
        0.,
        0.,
        0.,
        0.,
        1.092548430592079070,
        0.,
        #DZ2
        -0.315391565252520002,#*1.4142135623730951,  # extra sqrt(2),
        0.,
        0.,
        -0.315391565252520002,#*1.4142135623730951,  # extra sqrt(2),
        0.,
        0.630783130505040012,#*1.4142135623730951,  # extra sqrt(2),
        #DXZ
        0.,
        0.,
        1.092548430592079070,
        0.,  
        0., 
        0.,
        #DY2
#        0.546274215296039535*1.4142135623730951,  # extra sqrt(2),
#        1.092548430592079070, 
        0.546274215296039535,
#        1.092548430592079070/(2.**0.5), 
        0.,
        0.,
#        -1.092548430592079070,
#        -1.092548430592079070/(2.**0.5),
        -0.546274215296039535,
#        -0.546274215296039535*1.4142135623730951,  # extra sqrt(2),
        0.,
        0.,
        #F-3
        0.,
        1.770130769779930531,
        0.,
        0.,
        0.,
        0.,
        -0.590043589926643510,
        0.,
        0.,
        0.,
        #F-2
        0.,
        0.,
        0.,
        0.,
        2.890611442640554055,
        0.,
        0.,
        0.,
        0.,
        0.,
        #F-1
        0.,
        -0.457045799464465739,
        0.,
        0.,
        0.,
        0.,
        -0.457045799464465739,
        0.,
        1.828183197857862944,
        0.,
        #F0
        0.,
        0.,
        -1.119528997770346170,
        0.,
        0.,
        0.,
        0.,
        -1.119528997770346170,
        0.,
        0.746352665180230782,
        #F+1
        -0.457045799464465739,
        0.,
        0.,
        -0.457045799464465739,
        0.,
        1.828183197857862944,
        0.,
        0.,
        0.,
        0.,
        #F+2
        0.,
        0.,
        1.445305721320277020,
        0.,
        0.,
        0.,
        0.,
        -1.445305721320277020,
        0.,
        0.,
        #F+3
        0.590043589926643510,
        0.,
        0.,
        -1.770130769779930530,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        ]











