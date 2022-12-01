#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <vector>
#include <list>
#include <time.h>
#include <cstdlib>


    

template<class Kernel>
class Simulated_Annealing {
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private:   
        //internal helper functions
        double calculate_energy(double Area,double Hull_Area,std::string Criteria,int p_Size);
        bool validity_check(Triangle_2 t1,Triangle_2 t2,Polygon_2 Polygon);
        std::vector<std::vector<Point_2>> point_subsets(std::vector<Point_2>);
        void merge_polygons(std::vector<Polygon_2>);

        double local_step(Polygon_2& Polygon);
        double global_step(Polygon_2& Polygon);
        double sub_division(Polygon_2& Polygon);
    public:
        double solve(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,std::string Step_Choice,int Iterations);
        double solve(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,std::string Step_Choice,int Iterations,int Attempts);
        double expanded_solve(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,int Iterations,int Attempts=1);
};

#include "simulated_annealing.cpp"

#endif