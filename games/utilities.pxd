# distutils: sources = cpp/utilities.cpp
from libcpp.vector cimport vector
from CGames cimport MeanGame

cdef extern from "cpp/utilities.h":
    cdef struct LBF:
        vector[int] labeled_field;
        vector[int] cluster_sizes;

    ctypedef LBF LabeledField;

    vector[int] py_n_m_distribution(MeanGame*);
    LabeledField** clustering(const vector[int]&, int, int);
