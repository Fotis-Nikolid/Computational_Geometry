#ifndef SIMULATED_ANNEALING_CC
#define SIMULATED_ANNEALING_CC
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <list>
#include "simulated_annealing.h"
#include "incremental.h"
#include "hull.h"
#include <time.h>
#include <climits>


template<class Kernel> 
static bool comp_func(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2) {//used for topological sorting by x-axis in increasing order
    return (p1.x()<p2.x());
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
bool Simulated_Annealing<Kernel>::validity_check(Point_2 prev,Point_2 point1,Point_2 point2,Point_2 next,std::vector<Point_2> Points) {
    Triangle_2 t1,t2;
    t1=Triangle_2(prev,point1,point2);
    t2=Triangle_2(point1,point2,next);
    const auto cross_point=intersection(Segment_2(prev,point1),Segment_2(point2,next));//check if two new lines itersect with eachother
    if (cross_point) {
        return false;
    }
    for (auto point: Points) {
        //This is performed by checking if any point exists inside one of the two new triangles being formed by the new ordering of vertices
        //Has been discussed in class and shown to be a correct approach, which saves us the added complexity of doing a visibility check
        if (t1.bounded_side(point)!=CGAL::ON_UNBOUNDED_SIDE) {//if not outside of triangle, thus inside or on it's boundary, then polygon is not valid
            if (!(point==t1[0] || point==t1[1] || point==t1[2])) {
                return false;
            }
        }
        if (t2.bounded_side(point)!=CGAL::ON_UNBOUNDED_SIDE) {
            if (!(point==t2[0] || point==t2[1] || point==t2[2])) { 
                return false;
            }
        }
    }
    return true;
}

//Breaks the set of points into multiple subsets for which smaller polygon's will be created using global steps and then joined to then perform local steps
template<class Kernel>
std::vector<std::vector<CGAL::Point_2<Kernel>>> Simulated_Annealing<Kernel>::point_subsets(std::vector<Point_2> Points) {
    std::sort(Points.begin(), Points.end(), comp_func<Kernel>);
    int m=100;//average number of points per subset
    int k=std::ceil((Points.size()-1)/(m-1));//desirable number of subsets created, depending on the orientation of the points this is not always possible
    std::vector<std::vector<Point_2>> subsets;
    bool more_points=false;
    std::vector<Point_2> subset;
    for (int i=0;i<Points.size()-2;i++) {
        Point_2 p1,p2,p3;
        p1=Points.at(i);
        p2=Points.at(i+1);
        p3=Points.at(i+2);
        if (p1.y()<p2.y() && p3.y()<p2.y()) {//if next 3 consecutive points form the pattern needed for breaking into subsets
            if (i>((subsets.size()+1)*m*0.9) && (Points.size()-i)>(m*0.5)) {//only break if around the point expected of a subset break or afterwards
                subset.push_back(p1);
                subset.push_back(p2);//add the connecting point of the two subsets
                subsets.push_back(subset);
                subset.clear();

                subset.push_back(p2);
                subset.push_back(p3);
                i+=2;//skip two next points                
            }
            else {//if criteria for breaking is not met, simple add the point to current subset
                subset.push_back(p1);
            }
        }
        else {//if criteria for breaking is not met, simple add the point to current subset
            subset.push_back(p1);
        }
    }
    subset.push_back(Points.at(Points.size()-1));
    subset.push_back(Points.at(Points.size()-2));
    subsets.push_back(subset);  
    return subsets;
    
}
template<class Kernel>
CGAL::Polygon_2<Kernel> Simulated_Annealing<Kernel>::merge_polygons(std::vector<Polygon_2> polygons) {
    Polygon_2 merge_polygon=polygons.at(0);
    
    for (int i=1;i<polygons.size();i++) {
        Polygon_2 poly=polygons.at(i);
        Polygon_2 m=merge_polygon;

        typename Polygon_2::Vertices::iterator join_iter; //iterator to point that both polygons connect
        typename Polygon_2::Vertices::iterator prev_iter; //iterator to previous topological point to join point(it is downwards to the left)
        typename Polygon_2::Vertices::iterator next_iter; //iterator to next topological point to join point(it is downwards to the right)
        bool found_join=false;
        for (join_iter=merge_polygon.vertices_begin();join_iter<merge_polygon.vertices_end();join_iter++) {//find iterator of the point of connection between the two polygons
            for (auto poly_it=poly.vertices_begin();poly_it<poly.vertices_end();poly_it++) {
                if (*join_iter==*poly_it) {//find common point iterator between two polygons
                    next_iter=std::next(poly_it);//get next point to join point(belongs to the 2nd polygon)
                    if (next_iter==poly.vertices_end()) {//in case we loop around to the polygon's start
                        next_iter=poly.vertices_begin();
                    }
                    if (next_iter->y()>=poly_it->y()) {
                        if (poly_it==poly.vertices_begin()) {
                            next_iter=std::prev(poly.vertices_end());
                        }   
                        else {
                            next_iter=std::prev(poly_it);
                        }
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
        if (join_iter==merge_polygon.vertices_begin()) {
            prev_iter=std::prev(merge_polygon.vertices_end());
        }   
        else {
            prev_iter=std::prev(join_iter);
        }
        if (prev_iter->y()>=join_iter->y()) {
            prev_iter=std::next(join_iter);
            if (prev_iter==merge_polygon.vertices_end()) {//in case we loop around to the polygon's start
                prev_iter=merge_polygon.vertices_begin();
            }
        }
        
        //insert all points of the 2nd polygon after the previous point iterator, all the way to the join point(thus merging the two polygons) 
        int insert_prev=join_iter-merge_polygon.begin();
        auto it=next_iter;
        Point_2 join_point=*join_iter;
        std::vector<Point_2> v_points1,v_points2;
        bool next_vector=false;
        while (*it!=join_point) {//iterate from the next point iterator of the 2nd polygon all the way circling back until we arrive at the join point, inserting each polygons point along the way 
            if (!next_vector) {
                v_points1.push_back(*it);
            }
            else {
                v_points2.push_back(*it);
            }
            it=std::next(it);
            if (it==poly.vertices_end()) {
                next_vector=true;
                it=poly.vertices_begin();
            }
        }
        merge_polygon.insert(join_iter,v_points1.begin(),v_points1.end());
        insert_prev+=v_points1.size();
        if (v_points2.size()!=0) {
            merge_polygon.insert(merge_polygon.vertices_begin()+insert_prev,v_points2.begin(),v_points2.end());
        }

    }
    return merge_polygon;
}

//find 2 consecutive points and swap their positions
template<class Kernel>
bool Simulated_Annealing<Kernel>::local_step(Polygon_2& Polygon) {
    std::vector<int> choices;
    for (int i=0;i<Polygon.size();i++) {
        choices.push_back(i);
    } 
    while (choices.size()!=0) {
        Polygon_2 temp_Polygon(Polygon);
        int r_index=rand()%(choices.size());//pick a random point of the polygon, to begin the swap in the 2 following points
        int random_pick=choices[r_index];
        choices[r_index]=choices[choices.size()-1];
        choices.pop_back();
        typename Polygon_2::Vertices::iterator prev=temp_Polygon.vertices_begin();
        typename Polygon_2::Vertices::iterator swap1,swap2,next;
        Point_2 swap_t;
        
        prev=std::next(prev,random_pick);//point which is previous to the two consecutive points being swapped
        
        typename Polygon_2::Vertices::iterator iter=prev;
        for (int i=0;i<3;i++) {//get the iterators of the 3 following points
            iter=std::next(iter);
            if (iter==temp_Polygon.vertices_end()) {//given that some of the points might be after the end of the polygon(meaning we circle back), restart the iterator
                iter=temp_Polygon.vertices_begin();
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
        CGAL::Fuzzy_iso_box<CGAL::Search_traits_2<Kernel>> exact_range(Point_2(lower_x,lower_y),Point_2(upper_x,upper_y));
        tree.search( std::back_inserter(points), exact_range);
        
        if (validity_check(*prev,*swap1,*swap2,*next,points)) {//check if there new polygon is simple
            Polygon=temp_Polygon;
            return true;//return new polygon area and thus exit the loop
        }
        //if resulting polygon is not simple, the while loop will restart
    }
    return false;
    
}

template <class Kernel>
bool visible_points(const CGAL::Polygon_2<Kernel>& Polygon, const CGAL::Point_2<Kernel> &p1, const CGAL::Point_2<Kernel> &p2)
{
    CGAL::Segment_2<Kernel> seg = CGAL::Segment_2<Kernel>(p1, p2);
    for (CGAL::Segment_2<Kernel> PolygonEdge : Polygon.edges())
    {   
        if (seg==PolygonEdge || CGAL::Segment_2<Kernel>(p2, p1)==PolygonEdge ) {
            continue;
        }
        auto intersect_p = intersection(PolygonEdge, seg);
        if (intersect_p)
        {
            if (CGAL::Point_2<Kernel> *p = boost::get<CGAL::Point_2<Kernel>>(&*intersect_p))
            {
                if (*p != p1 && *p != p2)
                    return false;
            }
            else
            {
                return false;
            }
        }
    }

    return true;
}

template<class Kernel>
bool Simulated_Annealing<Kernel>::global_step(Polygon_2& Polygon,Point_2* no_change1,Point_2* no_change2) {
    //local search with L=1
    int insert_prev,insert_next;
    int move_pos,move_prev,move_next;
    Polygon_2 t_Poly;
    while (true) {
        bool loop=true;
        while (loop) { 
            loop=false;
            t_Poly=Polygon_2(Polygon);
            insert_prev = std::rand()%t_Poly.vertices().size();
            move_pos = std::rand()%t_Poly.vertices().size();
            int fail=0;
            int fail_check=Polygon.vertices().size();
            while(abs(insert_prev-move_pos)<=2) {
                fail++;
                if (fail>fail_check) return false;
                move_pos = std::rand()%t_Poly.vertices().size();
            }
        
        
            insert_next=insert_prev+1;
            move_next=move_pos+1;
            move_prev=(move_pos!=0)?(move_pos-1):(t_Poly.vertices().size()-1);
            if ((t_Poly.vertices_begin()+insert_next)==t_Poly.vertices_end()) {
                insert_next=0;
            }
            if ((t_Poly.vertices_begin()+move_next)==t_Poly.vertices_end()) {
                move_next=0;
            }

            if (no_change1!=NULL) {
                Point_2 p=*no_change1;
                if (p==*(t_Poly.vertices_begin()+move_pos)) loop=true;
                if (p==*(t_Poly.vertices_begin()+move_next)) loop=true;
                if (p==*(t_Poly.vertices_begin()+move_prev)) loop=true;
                if (p==*(t_Poly.vertices_begin()+insert_prev)) loop=true;
                if (p==*(t_Poly.vertices_begin()+insert_next)) loop=true;
            }
            if (no_change2!=NULL) {
                Point_2 p=*no_change2;
                if (p==*(t_Poly.vertices_begin()+move_pos)) loop=true;
                if (p==*(t_Poly.vertices_begin()+move_next)) loop=true;
                if (p==*(t_Poly.vertices_begin()+move_prev)) loop=true;
                if (p==*(t_Poly.vertices_begin()+insert_prev)) loop=true;
                if (p==*(t_Poly.vertices_begin()+insert_next)) loop=true;
                
            }
        }
        Point_2 moved_point=*(t_Poly.vertices_begin()+move_pos);

        t_Poly.erase((t_Poly.vertices_begin()+move_pos));
        if (move_pos<move_next) {
            move_next--;
        }
        if (move_pos<move_prev) {
            move_prev--;
        }
        if (move_pos<insert_prev) {
            insert_prev--;
        }
        if (move_pos<insert_next) {
            insert_next--;
        }
        
        Point_2 in_next=*(t_Poly.vertices_begin()+insert_next);
        Point_2 in_prev=*(t_Poly.vertices_begin()+insert_prev);
        Point_2 mv_next=*(t_Poly.vertices_begin()+move_next);
        Point_2 mv_prev=*(t_Poly.vertices_begin()+move_prev);
        
        
        t_Poly.insert((t_Poly.vertices_begin()+insert_next),moved_point);
        

        bool check1=visible_points<Kernel>(t_Poly,in_next,moved_point);
        bool check2=visible_points<Kernel>(t_Poly,in_prev,moved_point);
        bool check3=visible_points<Kernel>(t_Poly,mv_prev,mv_next);
        if (check1 && check2 && check3) {
            Polygon=t_Poly;
            break;
        }
  
    }   

    return true;
}

template<class Kernel>
bool Simulated_Annealing<Kernel>::sub_division(Polygon_2& Polygon,std::vector<Point_2> Points,std::string Criteria,int Iterations,double& initial_area) {
    srand((unsigned)time(NULL));
    std::vector<std::vector<Point_2>> Subsets;
    Subsets=point_subsets(Points);//break points into multiple subsets
    
    std::vector<Polygon_2> polygons;
    
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
        double init;
        Point_2 p1,p2;
        p1=subset[0];//first point in subset
        p2=subset[subset.size()-1];//last point in subset
        //use altered version of Hull algorithm which will make sure not to break any of the two lines needed for the criteria to work
        bool good=false;
        if (subset_n==0) {
            hull.solve(poly,points,crit,NULL,&p2);//only keep the last edge(first does not matter in the 1st subset)
            std::ofstream outfile;
            outfile.open("step"+std::to_string(subset_n)+".txt");
            solve(poly,Points,Criteria,"global",Iterations,init,NULL,&p2);
        }   
        else if (subset_n==(Subsets.size()-1)) {
            hull.solve(poly,points,crit,&p1,NULL);//only keep the first edge(last does not matter in the last subset)

            solve(poly,Points,Criteria,"global",Iterations,init,&p1,NULL);
        }
        else {
            hull.solve(poly,points,crit,&p1,&p2);//keep both edges
            solve(poly,Points,Criteria,"global",Iterations,init,&p1,&p2);
            
        }
        subset_n++;
        polygons.push_back(poly);
    } 

    //merge the polygons
    
    Polygon=merge_polygons(polygons);
    
    return abs(Polygon.area());//apply local step to the merged Polygon
    
}

//perform one attempt with local or global with either min max or random starting polygon
template<class Kernel>
double Simulated_Annealing<Kernel>::solve(Polygon_2& Polygon,std::vector<CGAL::Point_2<Kernel>> Points,std::string Criteria,std::string Step_Choice,int Iterations,double& initial_area,Point_2* no_change1,Point_2* no_change2) {
    srand((unsigned)time(NULL));
    if (!tree_exists) {
        for(auto point:Points) {
            tree.insert(point);//initialize a kd_tree from the set of points if it has not already been initialized
        }        
        tree_exists=true;
    }

    Polygon_2 convex_hull;
    CGAL::convex_hull_2(Points.begin(),Points.end(),std::back_inserter(convex_hull));
    double hull_area=abs(convex_hull.area());//needed to calculate state energy

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
    if (Polygon.size()==0) {
        Polygon=inc.getPolygon();
    }
    int size=Polygon.vertices().size();
    double area=inc.getPolygonArea();
    initial_area=area;
    
    double energy=calculate_energy(area,hull_area,Criteria,size);//calculate initial energy state of starting polygon
    
    Polygon_2 t_Polygon,best_Polygon(Polygon);
    double n_area,n_energy,p_area=area;

    double Temperature=1;
    bool success;
    bool no_change=true;
    int BadStepsCount=0;
    int BadStepsThreshold=30;
    int RetryCount=0;
    int RetryThreshold=2;
    int count=0;
    while (Temperature>0) {
        while (true) {
            t_Polygon=Polygon;
            if (Step_Choice=="local") {
                success=local_step(t_Polygon);
                n_area=abs(t_Polygon.area());
            }
            else if (Step_Choice=="global") {
                success=global_step(t_Polygon,no_change1,no_change2);
                count++;
                n_area=abs(t_Polygon.area());
            }
            if (!success) {//there exists no point combination capable of swapping
                std::cout<<"Deadlock, exiting..."<<std::endl;
                Temperature=-1;
                break;
            }
            n_energy=calculate_energy(n_area,hull_area,Criteria,size);//find energy of altered polygon
            
            double EnergyDifference=n_energy-energy;//energy difference
            if (EnergyDifference<0) {//if change is positive
                no_change=false;
                BadStepsCount++;
                if ((Criteria=="max" && n_area>p_area) || (Criteria=="min" && n_area<p_area)) {
                    BadStepsCount=0;
                    RetryCount=0;
                    p_area=n_area;
                    best_Polygon=t_Polygon;
                }
                break;//keep solution
            }
            else {//else if change is negative, only apply it by performing the Metropolis criteria
                BadStepsCount+=3;
                if (BadStepsCount>BadStepsThreshold) {//if we have taken too many bad actions without finding a better result, backtrack to best Polygon found and retry
                    t_Polygon=best_Polygon;
                    RetryCount++;
                    BadStepsCount=0;
                    //std::cout<<"Number of Bad Steps exceeded, Resetting..."<<std::endl;
                    break;
                }
                if (RetryCount>RetryThreshold) {
                    std::cout<<"Maximum Number of Resets reached, Stopping..."<<std::endl;
                    Temperature=-1;//if number of backtracks to the same best solution is exceeds threshold, stop loop and return best solution
                    break;
                }
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
    if (no_change) {
        std::cout<<"No optimizations found"<<std::endl;
    }
    Polygon=best_Polygon;
    area=p_area;
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
        Polygon_2 poly;
        //the starting polygon should be initialized as optimal area, while the rest should be initialized randomly
        area=solve(poly,Points,Criteria,Step_Choice,Iterations,t_area);
        if (Criteria=="max") {
            if (area>p_area) {
                p_area=area;
                initial_area=t_area;
                best_poly=poly;
            }
        }
        else {
            if (area<p_area) {
                p_area=area;
                initial_area=t_area;
                best_poly=poly;
            }
        }
    }
    Polygon=best_poly;
    return p_area;
}
#endif