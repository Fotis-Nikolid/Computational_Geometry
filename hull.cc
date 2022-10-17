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
    public:
        //internal helper functions
        bool overlaps_point(Triangle_2 triangle,std::vector<Point_2> Points) {
            typename std::vector<Point_2>::iterator point;
            
            for (point=Points.begin();point<Points.end();point++) {
                if (triangle.bounded_side(*point)!=CGAL::ON_UNBOUNDED_SIDE) {
                    return true;
                }
            }
            return false;
        }
        bool is_visible(Segment_2 Edge,Point_2 n_point,Polygon_2 Polygon) {
            Segment_2 n_edge1=Segment_2(Edge[0],n_point);
            Segment_2 n_edge2=Segment_2(Edge[1],n_point);
            for(const Segment_2& edge: Polygon.edges()){
                if (Edge==edge) {
                    continue;
                }
                if (intersection(n_edge1,edge)||intersection(n_edge2,edge)) {
                    return false;
                }
                
            }
            return true;
        }
        double Edge_Selection(Polygon_2 Polygon,Point_2 n_point,std::list<Point_2> remaining_points,std::string criteria) {
            Segment_2 selection;
            typename Polygon_2::Vertices::iterator it=Polygon.begin();
            int index,selection_index;
            switch(criteria) {
                case "1"://random choice of edge
                    index=0;
                    for(const Segment_2 edge : Polygon.edges()){
                        if (this.is_visible(edge,n_point,Polygon)) {
                            Triangle_2 triangle=Triangle_2(edge[0],edge[1],n_point);
                            if (this.overlaps_point(triangle,remaining_points)) {
                                index++;
                                continue;
                            }
                            
                            Polygon.insert(it+index,n_point);
                            return triangle.area();
                        }
                        else {
                            index++;
                        }

                    }
                    break;
                case "2"://pick edge that minimizes Area
                    double min_area=Polygon.area();
                    selection_index=0;
                    index=0;
                    for(const Segment_2 edge : Polygon.edges()){
                        if (this.is_visible(edge,n_point,Polygon)) {
                            Triangle_2 triangle=Triangle_2(edge[0],edge[1],n_point);
                            if (this.overlaps_point(triangle,remaining_points)) {
                                index++;
                                continue;
                            }
                            
                            Polygon.insert(it+index,n_point);
                            return triangle.area();

                        }
                        else {
                            index++;
                        }

                    }
                    break;
                case "3"://pick edge that maximizes Area
                    double max_area=0;
                    selection_index=0;
                    index=0;
                    for(const Segment_2 edge : Polygon.edges()){
                        if (this.is_visible(edge,n_point,Polygon)) {
                            Triangle_2 triangle=Triangle_2(edge[0],edge[1],n_point);
                            if (this.overlaps_point(triangle,remaining_points)) {
                                index++;
                                continue;
                            }
                            
                            Polygon.insert(it+index,n_point);
                            return triangle.area();
                        }
                        else {
                            index++;
                        }

                    }
                    break;
                default:
                    break;
            }
        } 
    public:
        double solve(Polygon_2& Polygon,std::list<Point_2>& Points,std::string Criteria) {
            typename std::list<Point_2>::iterator it;
            Point_2 n_point;
            double Area;

            CGAL::convex_hull_2(Points.begin(),Points.end(),std::back_inserter(Polygon));
            Area=Polygon.area();
            while (Points.size()>0) {
                n_point=Points.pop_back();
                Area-=Edge_Selection(Polygon,n_point,Points,Criteria);
            }
            return Area;
        }
        

};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Point_2<K>   Point;
typedef CGAL::Segment_2<K> Segment;

int main() {
    Hull<K> alg;
    //std::cout<<alg.solve(polygon,points,"1")<<std::endl;

    std::list<Point> points;
    Polygon result;
    points.push_back(Point(0,0));
    points.push_back(Point(4,0));
    points.push_back(Point(1,1));
    points.push_back(Point(3,1));
    points.push_back(Point(4,4));
    points.push_back(Point(4,0));



    return 0;

    
}