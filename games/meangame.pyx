# distutils: language = c++
import numpy as np
cimport cython

from libcpp.vector cimport vector
from numpy cimport import_array, PyArray_SimpleNewFromData, NPY_UINT8, NPY_UINT32, NPY_INT32, npy_intp, NPY_DOUBLE

from CGames cimport AbstractSpatialGame, TriangularFieldGame

cdef class TwoFieldSpatialGame:
    cdef:
        AbstractSpatialGame *c_game;
        cdef int _L, percfrom, perctill

    def __cinit__(self, int L, double b1, double b2, double lam, double mu, int percfrom=-1, perctill=-1, double K = 0, int seed = 42):
        self.c_game = new AbstractSpatialGame(L, b1, b2, lam, mu, K, seed)
        self._L = L
        self.percfrom = percfrom;
        self.perctill = perctill;

    def __dealloc__(self):
        del self.c_game

    @property
    def L(self):
        """"Return linear size of the field"""
        return self.c_game.size()
    
    @property
    def K(self):
        """ Return the measure of stochastic uncertainties (noise) """
        return self.c_game.get_K()
    
    @K.setter
    def K(self, new_K):
        """ Set measure of stochastic uncertainties """
        self.c_game.set_K(new_K)

    @property
    def field(self):
        """Return the field as a numpy array."""
        cdef npy_intp dims[1]
        dims[0] = 2*self._L*self._L
        f = PyArray_SimpleNewFromData(1, dims, NPY_UINT8, self.c_game.get_field_pointer())
        field1 = f[:self._L*self._L].reshape((1, self._L, self._L))
        field2 = f[self._L*self._L:].reshape((1, self._L, self._L))
        return np.concatenate((field1, field2), axis=0)

    @field.setter
    def field(self, arr):
        """Set the game field from 3 dimensional array"""
        arr = np.asarray(arr)
        if len(arr.shape) != 3:
            raise ValueError("Expected a 3D array, got %s-d." % len(arr.shape))
        if arr.size != 2*self._L*self._L:
            raise ValueError(f"Size mismatch: expected {2*self._L*self._L}, got {arr.size}.")

        ar1 = arr[0].ravel()
        ar2 = arr[1].ravel()
        cdef vector[int] field1, field2;
        field1.resize(self._L*self._L)
        field2.resize(self._L*self._L)
        for j in range(self._L*self._L):
            field1[j] = ar1[j]
            field2[j] = ar2[j]

        self.c_game.set_field(field1, field2)
        field1.clear()
        field2.clear()

    @property
    def b(self):
        cdef vector[double] bs = self.c_game.get_bs()
        b1 = bs[0]
        b2 = bs[1]
        return (b1, b2)

    @b.setter
    def b(self, arr):
        arr = np.asarray(arr)
        if len(arr.shape) != 1:
            raise ValueError("Expected a 1D array, got %s-d." % len(arr.shape))
        if arr.size != 2:
            raise ValueError(f"Size mismatch: expected 2 values, got {arr.size}.")
        self.c_game.set_b(arr[0], arr[1])

    @property
    def koef(self):
        cdef vector[double] koef = self.c_game.get_koef()
        lam = koef[0]
        mu = koef[1]
        return (lam, mu)

    @koef.setter
    def koef(self, arr):
        arr = np.asarray(arr)
        if len(arr.shape) != 1:
            raise ValueError("Expected a 1D array, got %s-d." % len(arr.shape))
        if arr.size != 2:
            raise ValueError(f"Size mismatch: expected 2 values, got {arr.size}.")
        self.c_game.set_koef(arr[0], arr[1])

    @property
    def densities(self):
        cdef npy_intp dims[1]
        dims[0] = self.c_game.get_densities_size()
        densities = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, self.c_game.get_densities_pointer())
        den1 = densities[0::2].reshape(1, -1)
        den2 = densities[1::2].reshape(1, -1)
        return np.concatenate((den1, den2), axis=0)

    @property
    def persistence(self):
        cdef vector[double] per = self.c_game.get_persistences()
        per1 = per[0]
        per2 = per[1]
        return (per1, per2)

    def mn_distribution(self):
        cdef vector[int] result = self.c_game.mn_distribution()
        cdef int x, y

        nmdist = np.zeros((2, 9, 9), dtype="int")
        for off in range(2):
            for x in range(9):
                for y in range(9):
                    nmdist[off, y, x] = result[off * 81 + y * 9 + x]
        return nmdist

    def evolve(self, int num_steps = 1):
        self.c_game.evolve(num_steps, self.percfrom, self.perctill)

cdef class TwoTriangularFieldSpatialGame(TwoFieldSpatialGame):

    def __cinit__(self, int L, double b1, double b2, double lam, double mu, int percfrom=-1, perctill=-1, double K = 0, int seed = 42):
        self.c_game = new TriangularFieldGame(L, b1, b2, lam, mu, K, seed)
        self._L = L
        self.percfrom = percfrom;
        self.perctill = perctill;


#################################################################
#                                                               #
#                       Utilities                               #
#                                                               #
#################################################################

cdef long[:, :] collors = np.array(((255, 255, 0, 255),
                                    (0, 0, 0, 255),
                                    (255, 0, 0, 255),
                                    (0, 0, 255, 255),
                                    (0, 0, 0, 255),
                                    (0, 255, 0, 255)), dtype=int)

cdef double[:, :] colors_f = np.array(((1, 1, 0, 1),
                                    (0, 0, 0, 1),
                                    (1, 0, 0, 1),
                                    (0, 0, 1, 1),
                                    (0, 0, 0, 1),
                                    (0, 1, 0, 1)), dtype=float)

@cython.cdivision(True)
@cython.boundscheck(False)
def _make_rgb(long[:, :] field):
    cdef:
        int L = field.shape[0]
        long[:, :, :] new_field = np.zeros((L, L, 3), dtype=int)

    for i in range(L):
        for j in range(L):
            for t in range(3):
                new_field[i, j, t] = collors[field[i, j]+2][t]

    return np.asarray(new_field)


def color_field_change(oldfield, newfield):
    return _make_rgb(newfield + 2*(newfield-oldfield))

import_array()