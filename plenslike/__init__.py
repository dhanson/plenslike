import os
import numpy as np
import ctypes as ct

from _plenslike import *

datadir = os.path.dirname(__file__) + "/data/"
pll     = ct.CDLL( os.path.dirname(__file__) + "/_plenslike.so")

# qest
class qest(ct.Structure):
    _fields_ = [ ("ntrm", ct.c_int),
                 ("lmax", ct.c_int),
                 ("s12L", ct.POINTER( ct.POINTER(ct.c_int) )),
                 ("w12L", ct.POINTER( ct.POINTER( ct.POINTER(ct.c_double) ))) ]

# mono
class plenslike_dat_mono(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmax",          ct.c_int),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("bl_fid",        ct.POINTER(ct.c_double)),
                 ("fl",            ct.POINTER(ct.c_double)),
                 ("vl_inv",        ct.POINTER(ct.c_double)),
                 ("al_inv",        ct.POINTER(ct.c_double)) ]

pll.load_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono), ct.c_char_p]
pll.free_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono) ]
pll.calc_plenslike_mono.restype        = ct.c_double
pll.calc_plenslike_mono_renorm.restype = ct.c_double

class mono():
    def __init__(self, fname):
        print "plenslike:: loading mono likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_mono()
        pll.load_plenslike_dat_mono( ct.byref(self.dat), fname)

    def __del__(self):
        pll.free_plenslike_dat_mono(self.dat)

    def calc_like(self, clpp):
        assert( len(clpp) >= self.dat.lmax )
        return pll.calc_plenslike_mono( ct.byref(self.dat),
                                        clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_like_renorm(self, clpp, cltt, bl):
        assert( len(clpp) >= self.dat.lmax )
        assert( len(cltt) >= self.dat.lmax )
        assert( len(bl)   >= self.dat.lmax )

        return pll.calc_plenslike_mono_renorm( ct.byref(self.dat),
                                               clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               cltt.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               bl.ctypes.data_as(   ct.POINTER(ct.c_double) ) )

    def calc_bins_clpp(self, clpp):
        bins = np.zeros( self.dat.nbins )

        pll.fill_plenslike_mono_bins( ct.byref(self.dat),
                                      clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                      bins.ctypes.data_as( ct.POINTER(ct.c_double) ) )

        return bins


# quad

# quad_qest
class plenslike_dat_quad_qest(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmax",          ct.c_int),
                 ("lmaxt",         ct.c_int),
                 ("lmax1",         ct.c_int),
                 ("lmax2",         ct.c_int),
                 ("lmax3",         ct.c_int),
                 ("lmax4",         ct.c_int),
                 ("s4hat",         ct.c_double),
                 ("s4std",         ct.c_double),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("vl_inv",        ct.POINTER(ct.c_double)),
                 ("rl_inv",        ct.POINTER(ct.c_double)),
                 ("sl_fid",        ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("bl1n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl2n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl3n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl4n1_fid",     ct.POINTER(ct.c_double)),
                 ("qe12",          ct.POINTER(qest)),
                 ("qe34",          ct.POINTER(qest)) ]

pll.load_plenslike_dat_quad_qest.argtypes   = [ ct.POINTER(plenslike_dat_quad_qest), ct.c_char_p]
pll.free_plenslike_dat_quad_qest.argtypes   = [ ct.POINTER(plenslike_dat_quad_qest) ]

pll.fill_quad_qest_resp_pp_blfid.argtypes   = [ ct.c_int, ct.POINTER(ct.c_double), ct.POINTER(plenslike_dat_quad_qest) ]


class quad_qest():
    def __init__(self, fname):
        print "plenslike:: loading quad_qest likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_quad_qest()
        pll.load_plenslike_dat_quad_qest( ct.byref(self.dat), fname)

    def __del__(self):
        pll.free_plenslike_dat_quad_qest(self.dat)

    def calc_qc_resp_pp_blfid(self, lmax, cltt):
        ret = np.zeros(lmax+1)
        pll.fill_quad_qest_resp_pp_blfid( lmax, ret.ctypes.data_as(ct.POINTER(ct.c_double)),
                                          ct.byref(self.dat), cltt.ctypes.data_as(ct.POINTER(ct.c_double)) )
        return ret

    def calc_like(self, clpp):
        pass

    def calc_like_renorm(self, clpp, cltt, bl):
        pass

    def calc_bins_clpp(self, clpp):
        pass
