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
//based on the criteria provided, picks an edge from the polygon to break and add the new point 
template<class Kernel>
double Hull<Kernel>::Edge_Selection(Polygon_2& Polygon,std::list<Point_2>& remaining_points,char criteria) {
    typename Polygon_2::Vertices::iterator it=Polygon.begin();
    int index,selection;
    double t_area_loss,min_area_loss,max_area_loss;
    Segment_2 max_loss_edge,min_loss_edge;

    
    std::vector<Point_2> point_candidates;
    
    bool visible_exists;
    std::vector<int> point_to_edge;//keep the index of the edge(from which we can acquire the indexes of it's two vertices) as to know the points-edges association 
    int position=0;
    std::unordered_map<Segment_2,int> placements;
    min_area_loss=std::numeric_limits<double>::max();
    max_area_loss=0;
    
    for(auto it=Polygon.edges().begin();it<Polygon.edges().end();++it){//iterate edges and for each one find the closest point, if there is one
        Segment_2 edge=*it;
        bool to_update=false;
        auto pair=pairings.find(edge);
        if (pair!=pairings.end() && last_inserted!=pair->second) {
            Triangle_2 triangle(pair->first[0],pair->first[1],pair->second);
            t_area_loss=abs(triangle.area());
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
        visible_exists=false;
        double min_dist=std::numeric_limits<double>::max();
        Point_2 closest_point;
        bool failure=true;
        for (const Point_2 point: remaining_points) {            
            double dist=CGAL::squared_distance(edge,point);
            if (dist<min_dist) {
                failure=false;
                min_dist=dist;
                closest_point=point;
            }
        }
        if (!failure && is_visible(edge,closest_point,Polygon)) {
            pairings[edge]=closest_point;
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
        else {
            pairings.erase(edge);
        }
        position++;
    }
    if (pairings.size()==0) {
        std::cout<<"No visible points found"<<std::endl;
        std::cout<<"Out of :"<<Polygon.edges().size()<<" candidates"<<std::endl;
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
            random=rand()%(pairings.size());
            
            access=pairings.begin();

            pair=*(std::next(access,random));
            new_point=pair.second;
            broken_edge=pair.first;

            placement=*placements.find(broken_edge);
            index=placement.second+1;
            //new_point=point_candidates.at(random);

            //index=point_to_edge.at(random);
            //broken_edge=Polygon.edge(index);
            triangle=Triangle_2(broken_edge[0],broken_edge[1],new_point);
            
            t_area_loss=abs(triangle.area());
            Polygon.insert(it+index,new_point);
            
            last_inserted=new_point;
            remaining_points.remove(new_point);
            
            pairings.erase(broken_edge);            
            return t_area_loss;
        case '2'://pick edge that maximizes Area loss, thus minimizing total polygon Area
            pair=*pairings.find(max_loss_edge);
            new_point=pair.second;
            broken_edge=max_loss_edge;

            placement=*placements.find(broken_edge);
            index=placement.second+1;
            
            Polygon.insert(it+index,new_point);
            
            last_inserted=new_point;
            remaining_points.remove(new_point);
            
            pairings.erase(broken_edge);            
            return max_area_loss;
        case '3'://pick edge that minimizes Area loss, maximizing total polygon Area
            pair=*pairings.find(min_loss_edge);
            new_point=pair.second;
            broken_edge=min_loss_edge;

            placement=*placements.find(broken_edge);
            index=placement.second+1;
            
            Polygon.insert(it+index,new_point);
            
            last_inserted=new_point;
            remaining_points.remove(new_point);
            
            pairings.erase(broken_edge);            
            return min_area_loss;
    }
    return 0;
} 

//only function call that is exposed to user
//given a list of Points and a criteria for edge selection, creates the Polygon and returns it's area 
template<class Kernel>
double Hull<Kernel>::solve(Polygon_2& Polygon,std::list<Point_2> Points,char Criteria) {
    Point_2 n_point;
    double Area;
    std::ofstream file;
    srand(time(0));
    /*
    file.open("steps.txt");
    for (const Point_2 p:Points) {
        file<<p<<std::endl;
    }
    file<<"-"<<std::endl;
    */
    CGAL::convex_hull_2(Points.begin(),Points.end(),std::back_inserter(Polygon));//initialize polygon from the convex hull, so that we can break it's edges and assimilate points into the polygon
    for (const Point_2 vertex: Polygon.vertices()) {
        Points.remove(vertex);//remove any point part of the convex hull from the list of points(as it is already part of the polygon)
        last_inserted=vertex;
    }
    Area=Polygon.area();
    while (Points.size()>0) {//iterate over list of points
        //std::cout<<Polygon.vertices().size()<<std::endl;
        double loss=Edge_Selection(Polygon,Points,Criteria);//add point to polygon by breaking of the appropriate edge and then update the polygon's area based on the area lost from assimilated a point
        if (loss==0) {
            return 0;
        }
        Area-=loss;
        /*
        for (const Segment_2 edge: Polygon.edges()) {
            file<<edge<<std::endl;
        }
        file<<"-"<<std::endl;
        */
    }
    return Area;
}

#endif