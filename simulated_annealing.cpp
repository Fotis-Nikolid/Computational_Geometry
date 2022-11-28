#ifndef SIMULATED_ANNEALING_CC
#define SIMULATED_ANNEALING_CC
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <vector>
#include <list>
#include "simulated_annealing.h"
#include "incremental.h"
#include <time.h>
#include <climits>


template<class Kernel> 
bool comp_func(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2) {
    return (p1.x()>p2.x());
}

template<class Kernel>
double calculate_energy(double Area,double Hull_Area,std::string Criteria,int p_Size) {
    if (Criteria=="max") {
        return p_Size*(1 - Area/Hull_Area);
    }
    else {
        return p_Size*(Area/Hull_Area);
    }
}     

//Checks whether after a change in order of two consecutive points, if the polygon remains simple
template<class Kernel>
bool Simulated_Annealing::validity_check(CGAL::Triangle_2<Kernel> t1,CGAL::Triangle_2<Kernel> t2,std::vector<CGAL::Point_2<Kernel>> Points) {
    for (auto point: Points) {
        //This is performed by checking if any point exists inside one of the two new triangles being formed by the new ordering of vertices
        //Has been discussed in class and shown to be a correct approach, which saves us the added complexity of doing a visibility check
        if (t1.bounded_side(point)!=CGAL::ON_UNBOUNDED_SIDE) {//if not outside of triangle, thus inside or on it's boundary, then polygon is not valid
            return false;
        }
        if (t2.bounded_side(point)!=CGAL::ON_UNBOUNDED_SIDE) {
            return false;
        }
    }
    return true;
}

//breaks the set of points into multiple subsets for which smaller polygon's will be created using global steps and then joined to then perform local step
template<class Kernel>
std::vector<std::vector<CGAL::Point_2<Kernel>>> Simulated_Annealing::point_subsets(std::vector<CGAL::Point_2<Kernel>> Points) {
    std::sort(Points.begin(), Points.end(), comp_func<Kernel>);
    int m=100;
    int k=std::ceil((Points.size()-1)/(m-1));
}

template<class Kernel>
double Simulated_Annealing::local_step(CGAL::Polygon_2<Kernel>& Polygon) {
    while (true) {
        int rand=rand()/Polygon.vertices().size();
        typename Polygon_2::Vertices::iterator prev,swap1,swap2,next;
        Point_2 swap_t;
        
        prev=Polygon.vertices().start();
        prev=std::next(prev,rand);
        
        auto iter=prev;
        for (int i=0;i<3;i++) {
            iter++;
            if (iter==Polygon.vertices().end()) {
                iter=polygon.vertices().start();
            }
            switch (i) {
                case 0:
                    swap1=iter;
                    break;
                case 1:
                    swap2=iter;
                    break;
                case 2:
                    next=iter;
                    break;
            }
        }

        swap_t=*swap2;
        *swap2=*swap1;
        *swap1=swap_t;

        std::vector<Point_2> points;
        //get points inside the rectangle
        //points=...

        Triangle_2 t1,t2;
        t1=Triangle_2(*prev,*swap1,*swap2);
        t2=Triangle_2(*swap1,*swap2,*next);
        if (!validity_check(t1,t2,points)) {
            return Polygon.area();
        }
    }
}   

template<class Kernel>
double Simulated_Annealing::global_step(CGAL::Polygon_2<Kernel>& Polygon) {
    //local search with L=1
}

template<class Kernel>
void Simulated_Annealing::sub_division(CGAL::Polygon_2<Kernel>& Polygon)3 {

}

template<class Kernel>
double Simulated_Annealing::solve(CGAL::Polygon_2<Kernel>& Polygon,std::vector<CGAL::Point_2<Kernel>> Points,std::string Criteria,std::string Step_Choice,int Iterations) {
    srand((unsigned)time(NULL));
    
    Polygon_2 convex_hull;
    CGAL::convex_hull_2(Points.begin(),Points.end(),std::back_inserter(convex_hull));
    double hull_area=convex_hull.area();//needed to calculate state energy

    char crit;
    if (Criteria=="max") {
        crit='3';
    }
    else if (Criteria=="min") {
        crit='2';
    }
    else {
        crit='1';
    }
    
    Incremental<Kernel> inc(Points,"1a",crit);//create a polygon which maximimes or minimizes total area
    Polygon=inc.getPolygon();
     
    int size=Polygon.vertices().size();
    double area=inc.getPolygonArea();
    
    double energy=calculate_energy(area,hull_area,Criteria,size);
    
    double Temperature=1;
    while (Temperature>=0) {
        while (true) {
            double n_area,n_energy;
            Polygon_2 t_Polygon=Polygon;
            if (Step_Choice=="local") {
                n_area=local_step(t_Polygon);
            }
            else if (Step_Choice=="global") {
                n_area=global_step(t_Polygon);
            }
            n_energy=calculate_energy(n_area,hull_area,Criteria,size);
            
            EnergyDifference=n_energy-energy;
            if (EnergyDifference<0) {
                break;
            }
            else {
                double R=(double rand())/RAND_MAX
                double threshold=std::exp(-EnergyDifference/Temperature);
                if (R<=threshold) {
                    break;
                }
            }
        }
        Polygon=t_Polygon;
        energy=n_energy;
        area=n_area;

        Temperature=Temperature-1.0/Iterations;
    }
    return area;
}

template<class Kernel>
double Simulated_Annealing::solve(CGAL::Polygon_2<Kernel>& Polygon,std::vector<CGAL::Point_2<Kernel>> Points,std::string Criteria,std::string Step_Choice,int Iterations,int Attempts) {
    Polygon_2 best_poly;
    double p_area=(Criteria=="max"?0:std::numeric_limits<double>::max());
    for (int i=0;i<Attempts;i++) {
        double area=solve(Polygon,Points,'random',Step_Choice,Iterations);
        if (Criteria=="max" && area>p_area) {
            if (area>p_area) {
                p_area=area;
                best_poly=Polygon;
            }
        }
        else {
            if (area<p_area) {
                p_area=area;
                best_poly=Polygon;
            }
        }
    }
    Polygon=best_poly;
    return p_area;
}
#endif