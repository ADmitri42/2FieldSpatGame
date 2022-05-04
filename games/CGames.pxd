from libcpp.vector cimport vector

cdef extern from "cpp/spatgame.h":
    cdef cppclass AbstractSpatialGame:
        AbstractSpatialGame(int, double, double, double, double);
        void evolve(int, int, int);
        vector[double] get_densities();
        int get_densities_size();
        int size();
        vector[int] get_field();
        void set_b(double, double);
        void set_koef(double, double);
        void set_field(vector[int], vector[int]);
        vector[double] get_bs();
        vector[double] get_koef();
        vector[double] get_persistences();
        vector[int] mn_distribution();

        double get_persistence();
        char* get_field_pointer();
        double* get_densities_pointer()
