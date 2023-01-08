from libcpp.vector cimport vector

cdef extern from "cpp/spatgame.h":
    cdef cppclass AbstractSpatialGame:
        AbstractSpatialGame(int, double, double, double, double, double, int);
        void evolve(int, int, int);
        vector[double] get_densities();
        int get_densities_size();
        int size();
        vector[int] get_field();
        double get_K();
        void set_b(double, double);
        void set_koef(double, double);
        void set_field(vector[int], vector[int]);
        void set_K(double);
        vector[double] get_bs();
        vector[double] get_koef();
        vector[double] get_persistences();
        vector[int] mn_distribution();

        double get_persistence();
        char* get_field_pointer();
        double* get_densities_pointer()
    
cdef extern from "cpp/triangular_field_game.h":
    cdef cppclass TriangularFieldGame(AbstractSpatialGame):
        TriangularFieldGame(int, double, double, double, double, double, int);
