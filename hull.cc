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
        bool overlaps_point(Triangle_2 triangle,std::list<Point_2> Points) {
            typename std::list<Point_2>::iterator it;
            
            for (const Point_2 point: Points) {
                if (triangle.bounded_side(point)!=CGAL::ON_UNBOUNDED_SIDE) {
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
                const auto cross_point1=intersection(n_edge1,edge);
                if (cross_point1) {
                    const Point_2* p=boost::get<Point_2 >(&*cross_point1);
                    if (*p!=(Point_2)Edge[0]) {
                        //std::cout<<"Cross:"<<*p<<" edge"<<(Point_2)Edge[0]<<std::endl;
                        return false;
                    }
                }
                const auto cross_point2=intersection(n_edge2,edge);
                if (cross_point2) {
                    const Point_2* p=boost::get<Point_2 >(&*cross_point2);
                    if (*p!=(Point_2)Edge[1]) {
                        //std::cout<<"Cross:"<<*p<<" edge"<<(Point_2)Edge[1]<<std::endl;
                        return false;
                    }
                }
                
            }
            return true;
        }
        double Edge_Selection(Polygon_2& Polygon,Point_2 n_point,std::list<Point_2> remaining_points,char criteria) {
            typename Polygon_2::Vertices::iterator it=Polygon.begin();
            int index,selection;
            double t_area_loss,min_area_loss,max_area_loss;
            switch(criteria) {
                case '1'://random choice of edge
                    index=1;
                    for(const Segment_2 edge : Polygon.edges()){
                        
                        if (is_visible(edge,n_point,Polygon)) {
                            Triangle_2 triangle=Triangle_2(edge[0],edge[1],n_point);
                            if (overlaps_point(triangle,remaining_points)) {
                                index++;
                                continue;
                            }
                            Polygon.insert(it+index,n_point);
                            return abs(triangle.area());
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
                            max_area_loss=t_area_loss;
                            selection=index;
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
        double solve(Polygon_2& Polygon,std::list<Point_2>& Points,char Criteria) {
            Point_2 n_point;
            double Area;

            CGAL::convex_hull_2(Points.begin(),Points.end(),std::back_inserter(Polygon));
            for (const Point_2 vertex: Polygon.vertices()) {
                Points.remove(vertex);
            }
            Area=Polygon.area();
            while (Points.size()>0) {
                n_point=Points.back();
                Points.pop_back();
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