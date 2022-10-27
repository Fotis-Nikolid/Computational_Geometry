#include <algorithm>
#include <limits>
#include "incremental.h"

template<class Kernel>
Incremental<Kernel>::Incremental(const std::vector<Point_2> Points, std::string how_to_sort, std::string how_to_remove_edge)
{
    std::vector<Point_2> ps(Points);
    Real_Polygon();

    this->sort(ps, how_to_sort);
    this->Initialize(ps, how_to_sort[0]);

    while(!ps.empty())
    {
        Point_2 point = ps.back();
        ps.pop_back();

        Incremental<Kernel>::RedEdgesBoundaries limits = this->find_red_edges_boundaries_and_recreate_convex_hull(point);

        this->construct_new_polygon(limits, point);
    }
}

template<class Kernel>
float Incremental<Kernel>::getPolygonArea()
{
    return Real_Polygon.area();
}

template<class Kernel>
CGAL::Polygon_2<Kernel> Incremental<Kernel>::getPolygon()
{
    return Real_Polygon;
}

template<class Kernel>
bool comp_x_less(CGAL::Point_2<Kernel> p1, CGAL::Point_2<Kernel> p2)
{
    if(p1.x == p2.x)
    {
        return (p1.y < p2.y);
    }
    return (p1.x < p2.x);
}
template<class Kernel>
bool comp_x_more(CGAL::Point_2<Kernel> p1, CGAL::Point_2<Kernel> p2)
{
    if(p1.x == p2.x)
    {
        return (p1.y < p2.y);
    }
    return (p1.x > p2.x);
}
template<class Kernel>
bool comp_y_less(CGAL::Point_2<Kernel> p1, CGAL::Point_2<Kernel> p2)
{
    if(p1.y == p2.y)
    {
        return (p1.x < p2.x);
    }
    return (p1.y < p2.y);
}
template<class Kernel>
bool comp_y_more(CGAL::Point_2<Kernel> p1, CGAL::Point_2<Kernel> p2)
{
    if(p1.y == p2.y)
    {
        return (p1.x < p2.x);
    }
    return (p1.y > p2.y);
}

template<class Kernel>
void Incremental<Kernel>::Sort(std::vector<Point_2> Points, std::string how_to_sort)
{
    if(how_to_sort == "1a")
    {
        std::sort(Points.begin(), Points.end(), comp_x_less);
    }
    else if(how_to_sort == "1b")
    {
        std::sort(Points.begin(), Points.end(), comp_x_more);
    }
    else if(how_to_sort == "2a")
    {
        std::sort(Points.begin(), Points.end(), comp_y_less);
    }
    else if(how_to_sort == "2b")
    {
        std::sort(Points.begin(), Points.end(), comp_y_more);
    }
}

template<class Kernel>
bool equal_three_points(CGAL::Point_2<Kernel> A, CGAL::Point_2<Kernel> B, CGAL::Point_2<Kernel> C, char x_or_y)
{
    if(x_or_y == '1')
    {
        return (A.x == B.x) && (B.x == C.x);
    }
    else
    {
        return (A.y == B.y) && (B.y == C.y);
    }
}

template<class Kernel>
bool equal_two_points(CGAL::Point_2<Kernel> A, CGAL::Point_2<Kernel> B, char x_or_y)
{
    if(x_or_y == '1')
    {
        return (A.x == B.x);
    }
    else
    {
        return (A.y == B.y);
    }
}

template<class Kernel>
void Incremental<Kernel>::Initialize(std::vector<Point_2> ps, char x_or_y)
{
    Point_2 A(ps.back());
    ps.pop_back();
    Point_2 B(ps.back());
    ps.pop_back();
    Point_2 C(ps.back());
    ps.pop_back();

    Real_Polygon.push_back(A);
    Real_Polygon.push_back(B);
    Real_Polygon.push_back(C);

    if(equal_three_points(A, B, C))
    {
        Point_2 p;
        do
        {
            p = Point_2(ps.back());
            Real_Polygon.push_back(p);
            ps.pop_back();            
        }while(equal_two_points(A, p));
    }

    Convex_Hull_Polygon(Real_Polygon);
}

template<class Kernel>
typename Incremental<Kernel>::RedEdgesBoundaries Incremental<Kernel>::find_red_edges_boundaries_and_recreate_convex_hull(Point_2 point)
{
    Incremental<Kernel>::RedEdgesBoundaries vertices;
    int iter_to_insert = 0;
    int edge_counter = 0;
    for(Segment_2 Convex_Hull_Edge : Convex_Hull_Polygon.edges())
    {
        //check if the edge is red
        if(this->visible(Convex_Hull_Edge, point))
        {
            if(iter_to_insert == 0)
            {
                iter_to_insert = edge_counter + 1;
                vertices.first_vertex = Convex_Hull_Edge[0];
                vertices.second_vertex = Convex_Hull_Edge[1];
            }
            else
            {
                vertices.second_vertex = Convex_Hull_Edge[1];
                Convex_Hull_Polygon.erase(Convex_Hull_Edge.vertices_begin() + edge_counter);
            }
        }
        else if(iter_to_insert != 0)
        {
            break;
        }
        edge_counter++;
    }

    Convex_Hull_Polygon.insert(Convex_Hull_Polygon.vertices_begin() + iter_to_insert);

    return vertices;
}

bool cmp_area_less(double area1, double area2)
{
    return area1 < area2 ;
}
bool cmp_area_more(double area1, double area2)
{
    return area1 > area2 ;
}

template<class Kernel>
void Incremental<Kernel>::construct_new_polygon(Incremental<Kernel>::RedEdgesBoundaries red_limits, Point_2 new_point, char how_to_remove_edge)
{
    typename Polygon_2::Vertices::iterator lower_limit_iter = Real_Polygon.vertices_begin();
    typename Polygon_2::Vertices::iterator uper_limit_iter = Real_Polygon.vertices_begin();
    int vertices_counter = 0;

    for(Point_2 p : Real_Polygon.vertices())
    {
        if(lower_limit_iter != -1 && uper_limit_iter != -1) break;

        if(p == red_limits.first_vertex)
        {
            lower_limit_iter += vertices_counter;

            if(uper_limit_iter != Real_Polygon.vertices_begin()) break;
        }
        else if(p == red_limits.second_vertex)
        {
            if(vertices_counter == 0)
            {
                uper_limit_iter = Real_Polygon.vertices_end();
            }
            else
            {
                uper_limit_iter += (vertices_counter - 1);
                break;
            }   
        }
        vertices_counter++;
    }

    double area;
    int (*compF)(double, double);

    switch(how_to_remove_edge)
    {
        case '1':
            return;
            break;
        case '2':
            area = std::numeric_limits<double>::max();
            compF = &cmp_area_more;
            break;
        case '3':
            area = 0.0;
            compF = &cmp_area_less;
            break;
    }

    Polygon_2 NewPolygon(Real_Polygon);

    for(typename Polygon_2::Vertices::iterator iter = lower_limit_iter ; iter <= uper_limit_iter ; iter++)
    {
        Polygon_2 temp_polygon(Real_Polygon);
        Segment_2 seg;

        if(iter == Real_Polygon.vertices_end())
        {
            seg = Segment_2(*iter, Real_Polygon.vertices_end());
        }
        else
        {
            seg = Segment_2(*iter, *(iter + 1));
        }

        if(this->visible(seg, new_point))
        {
            temp_polygon.inser(iter + 1, new_point);

            if(compF(area, temp_polygon.area()))
            {
                area = temp_polygon.area();
                NewPolygon = temp_polygon;
            }
        }
    }

    Real_Polygon = NewPolygon;
}

template<class Kernel>
bool Incremental<Kernel>::visible(Segment_2 seg, Point_2 new_point)
{
    for(Segment_2 PolygonEdge : Real_Polygon.edges())
    {
        if(PolygonEdge == seg)
        {
            continue;
        }

        for(int i = 0 ; i < 2 ; i++)
        {
            auto intersect_p = intersection(Segment_2(seg[i], new_point), PolygonEdge);

            if(intersect_p)
            {
                if(Point_2* p=boost::get<Point_2 >(&*intersect_p))
                {
                    if(*p!=(Point_2)seg[i]) return false;
                }
                else
                {
                    return false;
                }
            }
        }
    }

    return true;
}