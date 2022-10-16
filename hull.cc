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

        }
        double Edge_Selection(Polygon_2 Polygon,Point_2 n_point,std::list<Point_2> remaining_points,std::string criteria) {
            
        } 
    public:
        double solve(Polygon_2& Polygon,std::list<Point_2>& Points,std::string Criteria) {
            
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