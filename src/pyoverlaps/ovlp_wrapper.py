from ctypes import cdll, POINTER, c_long, c_double
import numpy as np
import os
import platform

# Include intwrap.so location on the  DL_LIBRARY_PATH 

def symmetrize(M):
    """Symmetryze an upper triangular matrix"""
    return(M + np.transpose(M) - np.diagflat(M.diagonal()))

def wrap_ovlp_matrix():
    pform = platform.system()
    if(pform == "Linux"):
        dll  = cdll.LoadLibrary("intwrap.so")
    elif(pform == "Windows"):
        dll = cdll.LoadLibrary("intwrap.dll")
    else:
        raise OSError("Platform " + pform + " currently not supported")
    func = dll.intOvlpCrt
    func.argtypes = [POINTER(c_long), POINTER(c_long),
                     POINTER(c_double), POINTER(c_double),
                     c_long, c_long]

    return(func)

def wrap_ovlp_matrix_mix():
    """Compute overlap matrix between two sets of 
    different basis functions"""
    pform = platform.system()
    if(pform == "Linux"):
        dll = cdll.LoadLibrary("intwrap.so")
    elif(pform == "Windows"):
        dll = cdll.LoadLibrary("intwrap.dll")
    else:
        raise OSError("Platform " + pform + " currently not supported")
    func = dll.intOvlpCrtMix
    func.argtypes = [POINTER(c_long), POINTER(c_long), POINTER(c_double),
                     POINTER(c_long), POINTER(c_long), POINTER(c_double),
                     POINTER(c_double), c_long, c_long, c_long, c_long]
    return(func)

def calc_cart_ovlp_matrix(bas1, atm1, env1, ncrt1,
                     bas2 = None, atm2 = None, env2 = None, ncrt2 = None):
    """Wrapper for the computation of the overlap matrix.
    All args are 1d flattened within the function"""
    nshl1 = len(bas1)
    Ncrt1 = ncrt1
    Ncrt2 = ncrt1
    bas1  = bas1.flatten().ctypes.data_as(POINTER(c_long))
    atm1  = atm1.flatten().ctypes.data_as(POINTER(c_long))
    env1  = env1.flatten().ctypes.data_as(POINTER(c_double))
    print("nshl1 = ", nshl1, "Ncrt1 = ", ncrt1)
    if(ncrt2 == None):
        __wrap_ovlp_matrix = wrap_ovlp_matrix()
    else:
        nshl2 = len(bas2)
        Ncrt2 = ncrt2
        bas2  = bas2.flatten().ctypes.data_as(POINTER(c_long))
        atm2  = atm2.flatten().ctypes.data_as(POINTER(c_long))
        env2  = env2.flatten().ctypes.data_as(POINTER(c_double))
        print("nshl2 = ", nshl2, "Ncrt2 = ", ncrt2)
        __wrap_ovlp_matrix = wrap_ovlp_matrix_mix()


    # Initialize Overlap matrix
    S = np.zeros((Ncrt1, Ncrt2)).flatten().ctypes.data_as(POINTER(c_double))

    # Call function
    if(ncrt2 == None):
        __wrap_ovlp_matrix(bas1, atm1, env1, S, c_long(Ncrt1), c_long(nshl1))
    else:
        __wrap_ovlp_matrix(bas1, atm1, env1, 
                               bas2, atm2, env2,
                               S, c_long(Ncrt1), c_long(nshl1),
                               c_long(Ncrt2), c_long(nshl2))

    S = np.reshape(S[:Ncrt1**2], (-1,Ncrt1)) 
    return(S)

