import os
import ctypes as ct

class plenslike_dat_basic(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmax",          ct.c_int),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("vl_inv",        ct.POINTER(ct.c_double)),
                 ("al_inv",        ct.POINTER(ct.c_double)),
                 ("fl",            ct.POINTER(ct.c_double)),
                 ("bl",            ct.POINTER(ct.c_double)) ]

plenslike_dat_basic_pointer = ct.POINTER(plenslike_dat_basic)
pll = ct.CDLL( os.path.dirname(__file__) + "/_plenslike.so")
pll.load_plenslike_dat_basic.restype = plenslike_dat_basic_pointer
pll.calc_plenslike_basic.restype     = ct.c_double

class basic():
    def __init__(self, tfname):
        self.dat = pll.load_plenslike_dat_basic(tfname)

    def __del__(self):
        pll.free_plenslike_dat_basic(self.dat)
