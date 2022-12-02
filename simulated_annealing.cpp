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
bool Simulated_Annealing<Kernel>::comp_func(const Point_2 p1, const Point_2 p2) {//used for topological sorting by x-axis in increasing order
    return (p1.x()>p2.x());
}

template<class Kernel>
double Simulated_Annealing<Kernel>::calculate_energy(double Area,double Hull_Area,std::string Criteria,int p_Size) {//"energy" of a polygon state following formula used in Simulated Annealing
    if (Criteria=="max") {
        return p_Size*(1 - Area/Hull_Area);
    }
    else {
        return p_Size*(Area/Hull_Area);
    }
}     

//Checks whether after a change in order of two consecutive points, if the polygon remains simple
template<class Kernel>
bool Simulated_Annealing<Kernel>::validity_check(Triangle_2 t1,Triangle_2 t2,std::vector<Point_2> Points) {
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
std::vector<std::vector<CGAL::Point_2<Kernel>>> Simulated_Annealing<Kernel>::point_subsets(std::vector<Point_2> Points) {
    std::sort(Points.begin(), Points.end(), comp_func);
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
template<class Kernel>
CGAL::Polygon_2<Kernel> Simulated_Annealing<Kernel>::merge_polygons(std::vector<Polygon_2> polygons) {
    Polygon_2 merge_polygon=polygons.at(0);
    for (int i=1;i<polygons.size();i++) {
        Polygon_2 poly=polygons.at(i);
        typename Polygon_2::Vertices::iterator join_iter; //iterator to point that both polygons connect
        typename Polygon_2::Vertices::iterator prev_iter; //iterator to previous topological point to join point(it is downwards to the left)
        typename Polygon_2::Vertices::iterator next_iter; //iterator to next topological point to join point(it is downwards to the right)
        bool found_join=false;
        for (join_iter=merge_polygon.vertices().start();join_iter<merge_polygon.vertices().end();join_iter++) {//find iterator of the point of connection between the two polygons
            for (auto poly_it=poly.vertices().start();poly_it<poly.vertices().end();poly_it++) {
                if (*join_iter==*poly_it) {//find common point iterator between two polygons
                    next_iter=poly_it;//get next point to join point(belongs to the 2nd polygon)
                    if (next_iter==poly.vertices().end()) {//in case we loop around to the polygon's start
                        next_iter=poly.vertices().start();
                    }
                    found_join=true;
                    break;
                    
                }
            }
            if (found_join) {
                break;
            }
        }
        //get previous point of join point(belongs to the 1st polygon), so that we connect the previous and the next points together(erasing their connection to the join point)
        if (join_iter==poly.vertices().start()) {
            prev_iter=std::prev(merge_polygon.vertices().end());
        }   
        else {
            prev_iter=std::prev(join_iter);
        }

        //insert all points of the 2nd polygon after the previous point iterator, all the way to the join point(thus merging the two polygons) 
        typename Polygon_2::Vertices::iterator insert_it=prev_iter;
        for (auto it=next_iter;*it!=*join_iter;it++) {//iterate from the next point iterator of the 2nd polygon all the way circling back until we arrive at the join point, inserting each polygons point along the way 
            merge_polygon.insert(insert_it,*it);//update the merged polygon
            insert_it=std::next(insert_it);//move the iterator along to continue inserting
        }
    }
}

//find 2 consecutive points and swap their positions
template<class Kernel>
double Simulated_Annealing<Kernel>::local_step(Polygon_2& Polygon) {
    while (true) {
        int random_pick=rand()/Polygon.vertices().size();//pick a random point of the polygon, to begin the swap in the 2 following points
        typename Polygon_2::Vertices::iterator prev,swap1,swap2,next;
        Point_2 swap_t;
        
        prev=Polygon.vertices().start();
        prev=std::next(prev,random_pick);//point which is previous to the two consecutive points being swapped
        
        auto iter=prev;
        for (int i=0;i<3;i++) {//get the iterators of the 3 following points
            iter++;
            if (iter==Polygon.vertices().end()) {//given that some of the points might be after the end of the polygon(meaning we circle back), restart the iterator
                iter=Polygon.vertices().start();
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
        std::vector<int> x_axis,y_axis;

        x_axis.push_back((*prev).x());
        x_axis.push_back((*swap1).x());
        x_axis.push_back((*swap2).x());
        x_axis.push_back((*next).x());

        y_axis.push_back((*prev).y());
        y_axis.push_back((*swap1).y());
        y_axis.push_back((*swap2).y());
        y_axis.push_back((*next).y());

        double lower_x=std::numeric_limits<double>::max();
        double lower_y=std::numeric_limits<double>::max();
        double upper_x=std::numeric_limits<double>::min();
        double upper_y=std::numeric_limits<double>::min();
    
        for (auto x:x_axis) {
            if (x<lower_x) lower_x=x;
            if (x>upper_x) upper_x=x;
        }
        for (auto y:y_axis) {
            if (y<lower_y) lower_y=y;
            if (y>upper_y) upper_y=y;
        }
        
        //in order to check if result creates a simple point, we must check whether any point exists inside of the two "triangles" being formed by the swap operation
        //however, as to minimize the amount of points this validity check should be performed
        //we use kd-trees to get only the points that exist inside the rectangle created by the four points
        kd_tree->find_points_inside_bounds(points, upper_x,upper_x, lower_x, lower_y);//upper x , y lower x , y (in this order)

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
double Simulated_Annealing<Kernel>::global_step(Polygon_2& Polygon) {
    //local search with L=1
}

template<class Kernel>
double Simulated_Annealing<Kernel>::sub_division(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,int Iterations) {
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
        
        Hull<Kernel> hull;
        std::list<Point_2> points(subset.begin(),subset.end());
        Segment_2 e1,e2;
        e1=Segment_2(subset[0],subset[1]);//first edge in subset
        e2=Segment_2(subset[subset.size()-1],subset[subset.size()-2]);//last edge in subset
        //use altered version of Hull algorithm which will make sure not to break any of the two lines needed for the criteria to work
        if (subset_n==0) {
            hull.solve(poly,points,crit,NULL,e2);//only keep the last edge(first does not matter in the 1st subset)
        }   
        else if (subset) {
            hull.solve(poly,points,crit,e1,NULL);//only keep the first edge(last does not matter in the last subset)
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
double Simulated_Annealing<Kernel>::solve(Polygon_2& Polygon,std::vector<CGAL::Point_2<Kernel>> Points,std::string Criteria,std::string Step_Choice,int Iterations,double& initial_area) {
    srand((unsigned)time(NULL));
    if (!tree_exists) {
        tree_exists=true;
        kd_tree=kdTree<Kernel>(Points);//initialize a kd_tree from the set of points if it has not already been initialized
    }

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
    initial_area=area;
    
    double energy=calculate_energy(area,hull_area,Criteria,size);//calculate initial energy state of starting polygon
    
    Polygon_2 t_Polygon;
    double n_area,n_energy;

    double Temperature=1;
    while (Temperature>=0) {
        while (true) {
            t_Polygon=Polygon;
            if (Step_Choice=="local") {
                n_area=local_step(t_Polygon);
            }
            else if (Step_Choice=="global") {
                n_area=global_step(t_Polygon);
            }
            n_energy=calculate_energy(n_area,hull_area,Criteria,size);//find energy of altered polygon
            
            double EnergyDifference=n_energy-energy;//energy difference
            if (EnergyDifference<0) {//if change is positive
                break;//keep solution
            }
            else {//else if change is negative, only apply it by performing the Metropolis criteria
                double R=((double) rand())/RAND_MAX;
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
double Simulated_Annealing<Kernel>::solve(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,std::string Step_Choice,int Iterations,int Attempts,double& initial_area) {
    Polygon_2 best_poly;
    double p_area=(Criteria=="max"?0:std::numeric_limits<double>::max());//initialize area depending on maximization or minimization
    for (int i=0;i<Attempts;i++) {
        double t_area;
        double area;
        if (i==0) {//only the first polygon should be initialized as optimal area, while the rest should be initialized randomly
            solve(Polygon,Points,Criteria,Step_Choice,Iterations,t_area);
        }
        else {
            solve(Polygon,Points,"random",Step_Choice,Iterations,t_area);
        }
        if (Criteria=="max" && area>p_area) {
            if (area>p_area) {
                p_area=area;
                initial_area=t_area;
                best_poly=Polygon;
            }
        }
        else {
            if (area<p_area) {
                p_area=area;
                initial_area=t_area;
                best_poly=Polygon;
            }
        }
    }
    Polygon=best_poly;
    return p_area;
}
#endif