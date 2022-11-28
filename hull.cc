#ifndef HULL_CC
#define HULL_CC
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <vector>
#include <list>
#include "hull.h"
#include <unordered_map>

    



//checks for visibility between the new point and the vertices of the edge we want to break
template<class Kernel>
bool Hull<Kernel>::is_visible(Segment_2 Edge,Point_2 n_point,Polygon_2 Polygon) {
    Segment_2 n_edge1=Segment_2(Edge[0],n_point);//create the two new edges that will be added to polygon if all goes well
    Segment_2 n_edge2=Segment_2(Edge[1],n_point);
    for(const Segment_2& edge: Polygon.edges()){//check if any of the two new edges intersects with any other edge
        if (Edge==edge) {
            continue;
        }
        const auto cross_point1=intersection(n_edge1,edge);
        if (cross_point1) {
            if (const Point_2* p=boost::get<Point_2 >(&*cross_point1)) {//intersection() can return either a Point or a Segment
                if (*p!=(Point_2)Edge[0]) {//if they intersect on a point other than their shared vertex, then it means that there is not visibility between the vertex and the new point
                    return false;
                }
            }
            else {//intersection() returned a semgent, meaning that the new point is exactly on top of one of polygon's edges, adding it is not allowed
                return false;
            }
        }
        //same as above
        const auto cross_point2=intersection(n_edge2,edge);
        if (cross_point2) {
            if (const Point_2* p=boost::get<Point_2 >(&*cross_point2)) {
                if (*p!=(Point_2)Edge[1]) {
                    return false;
                }
            }
            else {
                return false;
            }
        }
        
    }
    return true;
}
//finds every pair of edge-closest point, and based on the criteria, breaks a visible edge from the pair and add it's closest point to the polygon.
template<class Kernel>
double Hull<Kernel>::Edge_Selection(Polygon_2& Polygon,std::list<Point_2>& remaining_points,char criteria) {
    int index;
    double t_area_loss,min_area_loss,max_area_loss;
    Segment_2 max_loss_edge,min_loss_edge;
    
    std::unordered_map<Segment_2,int> placements;//for each edge we want to know it's position on the Polygon, so that we can know where to insert() the point
    
    //since we are iterating "all" the edges to find their closest point, we also calculate the min and max criteria so that we don't need to do a loop afterwards
    min_area_loss=std::numeric_limits<double>::max();
    max_area_loss=0;
    
    for(auto it=Polygon.edges().begin();it<Polygon.edges().end();++it){//iterate edges and for each one find the closest point, if there is one
        Segment_2 edge=*it;
        if (edge==first_edge || edge==last_edge) {//These two edges must never be broken
            continue;
        }
        auto pair=pairings.find(edge);//check if we have already found the closest point of the edge in a previous iteration
        
        //The Optimization!
        if (pair!=pairings.end() && last_inserted!=pair->second) {//if we already have found a point, and it is not the point that we just inserted to the polygon, then there is no point in wasting resources to find the same point again
            //we still want to find the area, for min max criteria and learn it's position in the new updated polygon 
            Triangle_2 triangle(pair->first[0],pair->first[1],pair->second);
            t_area_loss=abs(triangle.area());//to find the area that will be removed from the polygon, it is enough to check the area of the triangle "removed" from it
            if (t_area_loss<min_area_loss) {
                min_area_loss=t_area_loss;
                min_loss_edge=pair->first;
            }
            if (t_area_loss>max_area_loss) {
                max_area_loss=t_area_loss;
                max_loss_edge=pair->first;
            }
            placements.insert(std::pair<Segment_2,int>(edge,it-Polygon.edges().begin()));
            continue;
        }
        double min_dist=std::numeric_limits<double>::max();
        Point_2 closest_point;
        for (const Point_2 point: remaining_points) {//find the closest point for a new edge or an edge whose closest point was the last point inserted    
            double dist=CGAL::squared_distance(edge,point);
            if (dist<min_dist) {
                min_dist=dist;
                closest_point=point;
            }
        }
        if (is_visible(edge,closest_point,Polygon)) {//if the pair of an edge and it's closest point are visible, then "add" the pairing to the map
            pairings[edge]=closest_point;//inserts or updates the closest point for said edge in the map, meaning it will now be considered a candidate for possible breaking-inserting to the polygon
            placements.insert(std::pair<Segment_2,int>(edge,it-Polygon.edges().begin()));
            Triangle_2 triangle(edge[0],edge[1],closest_point);
            t_area_loss=abs(triangle.area());
            if (t_area_loss<min_area_loss) {
                min_area_loss=t_area_loss;
                min_loss_edge=edge;
            }
            if (t_area_loss>max_area_loss) {
                max_area_loss=t_area_loss;
                max_loss_edge=edge;
            }
        }
        else {//if the pair edge-closest point is not visible, then it is not a candidate for insertion
            pairings.erase(edge);//erase it completely from the map
        }
    }
    if (pairings.size()==0) {//if the are not visible points for any edge, then we reached a deadend, and as stated in the eclass thread, we print the error message and end the algorithm
        std::cout<<"No visible points found!!"<<std::endl;
        std::cout<<"Out of : "<<Polygon.edges().size()<<" edges"<<std::endl;
        std::cout<<"Num points left: "<<remaining_points.size()<<std::endl;
        std::cout<<"Criteria: "<<criteria<<std::endl;
        return 0;
    }
    Point_2 new_point;
    Triangle_2 triangle;
    Segment_2 broken_edge;
    int random;
    std::pair<Segment_2,Point_2> pair;
    std::pair<Segment_2,int> placement;
    typename std::unordered_map<Segment_2,Point_2>::iterator access; 
    switch(criteria) {
        case '1'://random choice of edges
            random=rand()%(pairings.size());//pick a random pair out of the map to break
            
            access=pairings.begin();

            pair=*(std::next(access,random));//get the random edge-point valid pair
            new_point=pair.second;
            broken_edge=pair.first;

            placement=*placements.find(broken_edge);//find the edge's position in the polygon
            index=placement.second+1;

            triangle=Triangle_2(broken_edge[0],broken_edge[1],new_point);
            
            t_area_loss=abs(triangle.area());//calculate area loss
            Polygon.insert(Polygon.begin()+index,new_point);//insert point to the correct index of the polygon, effectively breaking the edge and creating two new ones
            
            last_inserted=new_point;//used in informing the next iteration of the Edge_Selection to know for which edge-point pairs it is neccessary to recalculate the closest points
            remaining_points.remove(new_point);//remove from set of points left to insert
            
            pairings.erase(broken_edge);            
            return t_area_loss;
        case '2'://pick edge that maximizes Area loss, thus minimizing total polygon Area
            pair=*pairings.find(max_loss_edge);//get the pair that the edge that will achieve maximum area loss belongs
            new_point=pair.second;//closest point
            broken_edge=max_loss_edge;

            placement=*placements.find(broken_edge);
            index=placement.second+1;
            
            Polygon.insert(Polygon.begin()+index,new_point);
            
            last_inserted=new_point;
            remaining_points.remove(new_point);
            
            pairings.erase(broken_edge);            
            return max_area_loss;
        case '3'://pick edge that minimizes Area loss, maximizing total polygon Area
            pair=*pairings.find(min_loss_edge);//get the pair that the edge that will achieve minimum area loss belongs
            new_point=pair.second;
            broken_edge=min_loss_edge;

            placement=*placements.find(broken_edge);
            index=placement.second+1;
            
            Polygon.insert(Polygon.begin()+index,new_point);
            
            last_inserted=new_point;
            remaining_points.remove(new_point);
            
            pairings.erase(broken_edge);            
            return min_area_loss;
    }
    return 0;//never supposed to reach here...booo
} 

//only function call that is exposed
//given a list of Points and a criteria for edge selection, creates the Polygon and returns it's area 
template<class Kernel>
double Hull<Kernel>::solve(Polygon_2& Polygon,std::list<Point_2> Points,char Criteria,Edge_2 edge_to_keep1=NULL,Edge_2 edge_to_keep2=NULL) {
    Point_2 n_point;
    double Area;
    std::ofstream file;
    srand(time(0));
    
    first_edge=edge_to_keep1;
    last_edge=edge_to_keep2;

    #if 0  
    //used for testing purposes
        file.open("steps.txt");
        for (const Point_2 p:Points) {
            file<<p<<std::endl;
        }
        file<<"-"<<std::endl;
    #endif 

    CGAL::convex_hull_2(Points.begin(),Points.end(),std::back_inserter(Polygon));//initialize polygon from the convex hull, so that we can break it's edges and assimilate points into the polygon
    for (const Point_2 vertex: Polygon.vertices()) {
        Points.remove(vertex);//remove any point part of the convex hull from the list of points(as it is already part of the polygon)
        last_inserted=vertex;//for the Edge_Selection to work correctely without having to add a special if statement
    }
    Area=Polygon.area();
    while (Points.size()>0) {//iterate over list of points
        double loss=Edge_Selection(Polygon,Points,Criteria);//add point to polygon by breaking of the appropriate edge and then update the polygon's area based on the area lost from assimilated a point
        if (loss==0) {//an error has occured
            return 0;
        }
        Area-=loss;//reduce the polygon's area by the area of the triangle that was "cut off" from the polygon
        
        #if 0
            //used for testing purposes to show each step taken by the algorithm and make sure it works correctly
            //output file is then proccesed by a python script
            for (const Segment_2 edge: Polygon.edges()) {
                file<<edge<<std::endl;
            }
            file<<"-"<<std::endl;
        #endif
    }
    return Area;
}

#endif