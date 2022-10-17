#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <vector>
#include <list>

    

template<class Kernel>
class Hull {
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private: 
        float Area;

        //internal helper functions

        //given a Triangle, which essentialy is the combination of a polygon's edge and the point we want to add
        //check whether this would leave any of the remaining points outside of the new polygon
        bool overlaps_point(Triangle_2 triangle,std::list<Point_2> Points) {
            typename std::list<Point_2>::iterator it;
            
            for (const Point_2 point: Points) {
                //if a point is inside the triangle, then return true and pick another edge to break
                if (triangle.bounded_side(point)==CGAL::ON_BOUNDED_SIDE) {/
                    return true;
                }
            }
            return false;
        }
        //checks for visibility between the new point and the vertices of the edge we want to break
        bool is_visible(Segment_2 Edge,Point_2 n_point,Polygon_2 Polygon) {
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
        double Edge_Selection(Polygon_2& Polygon,Point_2 n_point,std::list<Point_2> remaining_points,char criteria) {
            typename Polygon_2::Vertices::iterator it=Polygon.begin();
            int index,selection;
            double t_area_loss,min_area_loss,max_area_loss;
            switch(criteria) {
                case '1'://random choice of edge
                    index=1;
                    for(const Segment_2 edge : Polygon.edges()){//iterate edges
                        if (is_visible(edge,n_point,Polygon)) {//check if new point has visibility to edge's vertices
                            Triangle_2 triangle=Triangle_2(edge[0],edge[1],n_point);
                            if (overlaps_point(triangle,remaining_points)) {//if adding point to polygon leaves other points outside, abort operation
                                index++;
                                continue;
                            }
                            Polygon.insert(it+index,n_point);//insert point between the two vertices of the edge(similar to breaking the edge and adding two new ones)
                            return abs(triangle.area());//return absolute area of the "triangle" we destroyed
                        }
                        else {
                            index++;
                        }

                    }
                    break;
                case '2'://pick edge that maximizes Area loss, thus minimizing total polygon Area
                    max_area_loss=0;
                    selection=1;
                    index=1;
                    for(const Segment_2 edge : Polygon.edges()){
                        
                        if (is_visible(edge,n_point,Polygon)) {
                            Triangle_2 triangle=Triangle_2(edge[0],edge[1],n_point);
                            t_area_loss=abs(triangle.area());
                            if (t_area_loss<=max_area_loss) {
                                index++;
                                continue;
                            }
                            if (overlaps_point(triangle,remaining_points)) {
                                index++;
                                continue;
                            }
                            max_area_loss=t_area_loss;//update maximum area loss we can achieve by breaking an edge
                            selection=index;//keep which edges causes the above maximum area loss
                            index++;
                            
                        }
                        else {
                            index++;
                        }

                    }
                    Polygon.insert(it+selection,n_point);
                    return max_area_loss;
                    
                    break;
                case '3'://pick edge that minimizes Area loss, maximizing total polygon Area
                    min_area_loss=Polygon.area();
                    selection=1;
                    index=1;
                    for(const Segment_2 edge : Polygon.edges()){
                        
                        if (is_visible(edge,n_point,Polygon)) {
                            Triangle_2 triangle=Triangle_2(edge[0],edge[1],n_point);
                            t_area_loss=abs(triangle.area());
                            if (t_area_loss>=min_area_loss) {
                                index++;
                                continue;
                            }
                            if (overlaps_point(triangle,remaining_points)) {
                                index++;
                                continue;
                            }
                            min_area_loss=t_area_loss;
                            selection=index;
                            index++;
                            

                        }
                        else {
                            index++;
                        }

                    }

                    Polygon.insert(it+selection,n_point);
                    return min_area_loss;
                    break;
                default:
                    break;
            }
        } 
    public:
        //only function call that is exposed to user
        //given a list of Points and a criteria for edge selection, creates the Polygon and returns it's area 
        double solve(Polygon_2& Polygon,std::list<Point_2>& Points,char Criteria) {
            Point_2 n_point;
            double Area;

            CGAL::convex_hull_2(Points.begin(),Points.end(),std::back_inserter(Polygon));//initialize polygon from the convex hull, so that we can break it's edges and assimilate points into the polygon
            for (const Point_2 vertex: Polygon.vertices()) {
                Points.remove(vertex);//remove any point part of the convex hull from the list of points(as it is already part of the polygon)
            }
            Area=Polygon.area();
            while (Points.size()>0) {//iterate over list of points
                n_point=Points.back();
                Points.pop_back();
                Area-=Edge_Selection(Polygon,n_point,Points,Criteria);//add point to polygon by breaking of the appropriate edge and then update the polygon's area based on the area lost from assimilated a point
                
            }
            return Area;
        }
        

};
//testing below
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Point_2<K>   Point;
typedef CGAL::Segment_2<K> Segment;

int main() {
    Hull<K> alg;

    std::list<Point> points;
    Polygon result;
    points.push_back(Point(0,0));
    points.push_back(Point(4,0));
    
    points.push_back(Point(4,4));
    points.push_back(Point(0,4));
    points.push_back(Point(1,1));
    points.push_back(Point(1,3));

    std::cout<<alg.solve(result,points,'3')<<std::endl;

    for(const Segment edge : result.edges()){
        std::cout<<edge<<std::endl;
    }

    return 0;

    
}