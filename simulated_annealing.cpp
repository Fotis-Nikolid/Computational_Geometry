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
#include "hull.h"
#include <time.h>
#include <climits>


template<class Kernel> 
bool comp_func(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2) {//used for topological sorting by x-axis in increasing order
    return (p1.x()>p2.x());
}

template<class Kernel>
double calculate_energy(double Area,double Hull_Area,std::string Criteria,int p_Size) {//"energy" of a polygon state following formula used in Simulated Annealing
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

//Breaks the set of points into multiple subsets for which smaller polygon's will be created using global steps and then joined to then perform local steps
template<class Kernel>
std::vector<std::vector<CGAL::Point_2<Kernel>>> Simulated_Annealing::point_subsets(std::vector<CGAL::Point_2<Kernel>> Points) {
    std::sort(Points.begin(), Points.end(), comp_func<Kernel>);
    int m=100;//average number of points per subset
    int k=std::ceil((Points.size()-1)/(m-1));//desirable number of subsets created, depending on the orientation of the points this is not always possible
    std::vector<std::vector<Point_2>> subsets;
    int subset_n=1;
    for (int i=0;i<Points.size()-2;i++) {
        Point_2 p1,p2,p3;
        p1=Points.at(i);
        p2=Points.at(i+1);
        p3=Points.at(i+2);
        if (p1.y()<p2.y() && p3.y()<p2.y()) {//if next 3 consecutive points form the pattern needed for breaking into subsets
            if (i>(subset_n*m*0.9)) {//only break if around the point expected of a subset break or afterwards
                subsets[subset_n].push_back(p1);
                subsets[subset_n].push_back(p2);//add the connecting point of the two subsets

                //subsets[subset_n+1].append(p2); //these two operations will occure naturally in the next two loops
                //subsets[subset_n+1].append(p3);
                subset_n++;//following iterations will add points  
            }
        }
        else {//if criteria for breaking is not met, simple add the point to current subset
            subsets[subset_n].push_back(p1);
        }
        
    }
    return subsets;
    
}

//find 2 consecutive points and swap their positions
template<class Kernel>
double Simulated_Annealing::local_step(CGAL::Polygon_2<Kernel>& Polygon) {
    while (true) {
        int rand=rand()/Polygon.vertices().size();//pick a random point of the polygon, to begin the swap in the 2 following points
        typename Polygon_2::Vertices::iterator prev,swap1,swap2,next;
        Point_2 swap_t;
        
        prev=Polygon.vertices().start();
        prev=std::next(prev,rand);//point which is previous to the two consecutive points being swapped
        
        auto iter=prev;
        for (int i=0;i<3;i++) {//get the iterators of the 3 following points
            iter++;
            if (iter==Polygon.vertices().end()) {//given that some of the points might be after the end of the polygon(meaning we circle back), restart the iterator
                iter=polygon.vertices().start();
            }
            switch (i) {
                case 0:
                    swap1=iter;//1st follwing point to be swapped 
                    break;
                case 1:
                    swap2=iter;//2nd following point to be swapped with 1st
                    break;
                case 2:
                    next=iter;//point after the 2 that will be swapped
                    break;
            }
        }
        //perform the swap
        swap_t=*swap2;
        *swap2=*swap1;
        *swap1=swap_t;

        std::vector<Point_2> points;
        
        //in order to check if result creates a simple point, we must check whether any point exists inside of the two "triangles" being formed by the swap operation
        //however, as to minimize the amount of points this validity check should be performed
        //we use kd-trees to get only the points that exist inside the rectangle created by the four points
        
        //points=...

        Triangle_2 t1,t2;
        t1=Triangle_2(*prev,*swap1,*swap2);
        t2=Triangle_2(*swap1,*swap2,*next);
        if (!validity_check(t1,t2,points)) {//check if there new polygon is simple
            return Polygon.area();//return new polygon area and thus exit the loop
        }
        //if resulting polygon is not simple, the while loop will restart
    }
}   

template<class Kernel>
double Simulated_Annealing::global_step(CGAL::Polygon_2<Kernel>& Polygon) {
    //local search with L=1
}

template<class Kernel>
void Simulated_Annealing::sub_division(CGAL::Polygon_2<Kernel>& Polygon,std::vector<CGALL::Point_2<Kernel>> Points,std::string Criteria,int Iterations) {
    std::vector<std::vector<Point_2>> Subsets;
    
    Subsets=point_subsets(Points);//break points into multiple subsets
    
    std::list<Polygon_2> polygons;
    
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
    
    int subset_n=0;
    for (auto subset:Subsets) {
        Polygon_2 poly;
        
        Hull hull;
        std::list<Point_2> points(subset.begin(),subset.end());
        Edge_2 e1,e2;
        e1=Edge_2(subset[0],subset[1]);//first edge in subset
        e2=Edge_2(subset[subset.size()-1],subset[subset.size()-2]);//last edge in subset
        //use altered version of Hull algorithm which will make sure not to break any of the two lines needed for the criteria to work
        if (subset_n==0) {
            hull.solve(poly,points,crit,NULL,e2);//only keep the last edge(first does not matter in the 1st subset)
        }   
        else if (subset) {
            hull.solve(poly,points,crit,e1,NULL)//only keep the first edge(last does not matter in the last subset)
        }
        else {
            hull.solve(poly,points,crit,e1,e2);//keep both edges
        }
        solve(poly,Criteria,"global",Iterations);//solve subset based on global 

        polygons.append(poly);
        
    } 
    //merge the polygons
    //apply local step 
    
}

//perform one attemp with local or global with either min max or random starting polygon
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
    
    Incremental<Kernel> inc(Points,"1a",crit);//create a polygon which maximimes/minimizes/randomizes total area
    Polygon=inc.getPolygon();
     
    int size=Polygon.vertices().size();
    double area=inc.getPolygonArea();
    
    double energy=calculate_energy(area,hull_area,Criteria,size);//calculate initial energy state of starting polygon
    
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
            n_energy=calculate_energy(n_area,hull_area,Criteria,size);//find energy of altered polygon
            
            EnergyDifference=n_energy-energy;//energy difference
            if (EnergyDifference<0) {//if change is positive
                break;//keep solution
            }
            else {//else if change is negative, only apply it by performing the Metropolis criteria
                double R=(double rand())/RAND_MAX
                double threshold=std::exp(-EnergyDifference/Temperature);//Metropolis formula
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

//perform mulitple attempts with randomized starting polygons and pick the best polygon out of all
template<class Kernel>
double Simulated_Annealing::solve(CGAL::Polygon_2<Kernel>& Polygon,std::vector<CGAL::Point_2<Kernel>> Points,std::string Criteria,std::string Step_Choice,int Iterations,int Attempts) {
    Polygon_2 best_poly;
    double p_area=(Criteria=="max"?0:std::numeric_limits<double>::max());//initialize area depending on maximization or minimization
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