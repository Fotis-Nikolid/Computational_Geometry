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
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>


    

template<class Kernel>
class Simulated_Annealing {
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private:   
        CGAL::Kd_tree<CGAL::Search_traits_2<Kernel>> tree;
        bool tree_exists=false;
        //internal helper functions
        double calculate_energy(double Area,double Hull_Area,std::string Criteria,int p_Size);
        bool validity_check(Point_2,Point_2,Point_2,Point_2,std::vector<Point_2> Points);
        std::vector<std::vector<Point_2>> point_subsets(std::vector<Point_2>);
        Polygon_2 merge_polygons(std::vector<Polygon_2>);

        bool local_step(Polygon_2& Polygon);
        bool global_step(Polygon_2& Polygon,Point_2*,Point_2*);
        
    public:
        double solve(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,std::string Step_Choice,int Iterations,double& Init_Area,Point_2* no_change1=NULL,Point_2* no_change2=NULL);
        double solve(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,std::string Step_Choice,int Iterations,int Attempts,double& Init_Area);
        bool sub_division(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,int Iterations,double& initial_area);
};

#include "simulated_annealing.cpp"

#endif