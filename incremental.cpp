#ifndef INCREMENTAL_CPP
#define INCREMENTAL_CPP
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "incremental.h"


template<class Kernel> Incremental<Kernel>::Incremental(const std::vector<Point_2> Points, const std::string how_to_sort, const char how_to_remove_edge)
{
    std::vector<Point_2> ps(Points);
    //sort the ps vector
    this->Sort(ps, how_to_sort);
    //create the triangle (if the have the same x or y then its not a triangle its a polygon)
    this->Initialize(ps, how_to_sort[0]);

    //keep removing point from the end of the vector
    while(!ps.empty())
    {
        Point_2 point = ps.back();
        ps.pop_back();
        //find the red edges and create the new convex hull polygon with the new point
        Incremental<Kernel>::RedEdgesBoundaries limits = this->find_red_edges_boundaries_and_recreate_convex_hull(point);
        //find a edge of the real polygon and add the new point
        this->construct_new_polygon(limits, point, how_to_remove_edge);
    }
}

template<class Kernel> float Incremental<Kernel>::getPolygonArea()
{
    return std::abs(Real_Polygon.area());
}

template<class Kernel> CGAL::Polygon_2<Kernel> Incremental<Kernel>::getPolygon()
{
    return Real_Polygon;
}

template<class Kernel> bool Incremental<Kernel>::comp_x_less(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2)
{
    if(p1.x() == p2.x())
    {
        return (p1.y() < p2.y());
    }
    return (p1.x() < p2.x());
}
template<class Kernel> bool Incremental<Kernel>::comp_x_more(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2)
{
    if(p1.x() == p2.x())
    {
        return (p1.y() < p2.y());
    }
    return (p1.x() > p2.x());
}
template<class Kernel> bool Incremental<Kernel>::comp_y_less(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2)
{
    if(p1.y() == p2.y())
    {
        return (p1.x() < p2.x());
    }
    return (p1.y() < p2.y());
}
template<class Kernel> bool Incremental<Kernel>::comp_y_more(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2)
{
    if(p1.y() == p2.y())
    {
        return (p1.x() < p2.x());
    }
    return (p1.y() > p2.y());
}

template<class Kernel> void Incremental<Kernel>::Sort(std::vector<Point_2>& Points, const std::string how_to_sort)
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

template<class Kernel> bool equal_three_points(const CGAL::Point_2<Kernel> A, const CGAL::Point_2<Kernel> B, const CGAL::Point_2<Kernel> C, const char x_or_y)
{
    if(x_or_y == '1')
    {
        return (A.x() == B.x()) && (B.x() == C.x());
    }
    else
    {
        return (A.y() == B.y()) && (B.y() == C.y());
    }
}

template<class Kernel> bool equal_two_points(const CGAL::Point_2<Kernel> A, const CGAL::Point_2<Kernel> B, const char x_or_y)
{
    if(x_or_y == '1')
    {
        return (A.x() == B.x());
    }
    else
    {
        return (A.y() == B.y());
    }
}

template<class Kernel> void Incremental<Kernel>::Initialize(std::vector<Point_2>& ps, const char x_or_y)
{
    //and the last 3 points of the vecotr in the trinagle and remove them from the vector
    Point_2 A(ps.back());
    ps.pop_back();
    Point_2 B(ps.back());
    ps.pop_back();
    Point_2 C(ps.back());
    ps.pop_back();
    
    Real_Polygon.push_back(A);
    Real_Polygon.push_back(B);
    Real_Polygon.push_back(C);

    //if the 3 points have the same x or y (depending if are sorting with x or y)
    //then keep adding points to the polygon until there are different
    if(equal_three_points(A, B, C, x_or_y))
    {
        Point_2 p;
        do
        {
            p = Point_2(ps.back());
            Real_Polygon.push_back(p);
            ps.pop_back();            
        }while(equal_two_points(A, p, x_or_y));
    }
    
    Convex_Hull_Polygon = Real_Polygon;
}

template<class Kernel> typename Incremental<Kernel>::RedEdgesBoundaries Incremental<Kernel>::find_red_edges_boundaries_and_recreate_convex_hull(const Point_2 point)
{
    Incremental<Kernel>::RedEdgesBoundaries vertices;
    int iter_to_insert = 0;
    int edge_counter = 0;

    std::vector<int> Vertices_to_remove;//inside vertices that we be removed to create the new convex_hull_polygon

    //the red edges are continous in the convex hull polygon
    //we find the the vertex the red edges start (we save it as vertices.first_vertex)
    //and the vertex that red edges end (we save it as vertices.second_vertex)
    for(Segment_2 Convex_Hull_Edge : Convex_Hull_Polygon.edges())
    {
        //check if the edge is red
        if(this->red_visible(Convex_Hull_Edge, point))
        {
            //if this is the first red edge we find
            if(iter_to_insert == 0)
            {
                iter_to_insert = edge_counter + 1;
                vertices.first_vertex = Convex_Hull_Edge[0];
                vertices.second_vertex = Convex_Hull_Edge[1];
            }
            //if its not the first red edge we find
            else
            {
                vertices.second_vertex = Convex_Hull_Edge[1];
                //we need to remove Convex_Hull_Edge[0] to create the new convex hull
                Vertices_to_remove.push_back(edge_counter);
            }
        }
        //if we have already found red edges and this edge is not red then there are no more red edges in the polygon
        else if(iter_to_insert != 0)
        {
            break;
        }
     
        edge_counter++;
    }

    //create the new convex hull
    Convex_Hull_Polygon.insert(Convex_Hull_Polygon.vertices_begin() + iter_to_insert, point);

    int erased_elm = 0;

    for(std::vector<int>::iterator iter = Vertices_to_remove.begin() ; iter < Vertices_to_remove.end() ; iter++)
    {
        Convex_Hull_Polygon.erase(Convex_Hull_Polygon.vertices_begin() + (*iter) - erased_elm + 1);
        erased_elm++;
    }

    return vertices;
}

template<class Kernel> void Incremental<Kernel>::construct_new_polygon(const Incremental<Kernel>::RedEdgesBoundaries red_limits, const Point_2 new_point, const char how_to_remove_edge)
{
    typename Polygon_2::Vertices::iterator lower_limit_iter = Real_Polygon.vertices_begin();
    typename Polygon_2::Vertices::iterator uper_limit_iter = Real_Polygon.vertices_begin();
    int vertices_counter = 0;

    //find the range of the candidate edges that can be broken and connected with the new vertex
    //we are searching for first_vertex and second_vertex of the red_limits
    //expect in one case , the first_vertex will be always found before the second_vertex
    //in the case that the last edge of the polygon is in that range then the second_vertex will be found in vertices_begin()
    for(Point_2 p : Real_Polygon.vertices())
    {
        //if we have found the range then there no point in looking up more vertices
        if(lower_limit_iter != Real_Polygon.vertices_begin() && uper_limit_iter != Real_Polygon.vertices_begin()) break;

        //if we found the first_vertex
        if(p == red_limits.first_vertex)
        {
            //save the iterator
            lower_limit_iter += vertices_counter;
            //if we are at the case that the second_vertex was found in potition 0 then there no reason to continue the loop
            if(uper_limit_iter != Real_Polygon.vertices_begin()) break;
        }
        //if we found the second_vertex
        else if(p == red_limits.second_vertex)
        {
            //if we found it in potition 0 then set the uper_limit_iter to the last vertex of the polygon
            if(vertices_counter == 0)
            {
                uper_limit_iter = Real_Polygon.vertices_end() - 1;
            }
            else
            {
                //else just set the uper_limit_iter to the vertex before of that we found
                uper_limit_iter += (vertices_counter - 1);
                //if we did not found the vertex at potiton 0 it means we have arleady found the first_vertex too so we stop the loop
                break;
            }   
        }
        vertices_counter++;
    }

    //in the case we are picking random edge to brake
    if(how_to_remove_edge == '1')
    {
        std::srand(std::time(nullptr));
        
        std::vector<int> Visible_Edges;//vector that saves the potition of the visible edges
        //from the candidate edges find the visible ones
        for(typename Polygon_2::Vertices::iterator iter = lower_limit_iter ; iter <= uper_limit_iter ; iter++)
        {
            Segment_2 seg;

            if(iter == Real_Polygon.vertices_end() - 1)
            {
                seg = Segment_2(*iter, *(Real_Polygon.vertices_begin()));
            }
            else
            {
                seg = Segment_2(*iter, *(iter + 1));
            }

            if(this->visible(seg, new_point))
            {
                Visible_Edges.push_back(iter - Real_Polygon.begin());
            }
        }
        //pick a random visible edge to break
        int random_pick = std::rand()%(Visible_Edges.size());
        //insert the two new edges in the random potition
        Real_Polygon.insert(Visible_Edges[random_pick] + 1 + Real_Polygon.begin(), new_point);
    }
    //if we try to find the min area polygon
    else if(how_to_remove_edge == '2')
    {
        double area = std::numeric_limits<double>::max();
        Polygon_2 NewPolygon(Real_Polygon);

        for(typename Polygon_2::Vertices::iterator iter = lower_limit_iter ; iter <= uper_limit_iter ; iter++)
        {
            Triangle_2 temp_triangle;
            Segment_2 seg;

            if(iter == Real_Polygon.vertices_end() - 1)
            {
                temp_triangle = Triangle_2(*iter, new_point, *(Real_Polygon.vertices_begin()));
                seg = Segment_2(*iter, *(Real_Polygon.vertices_begin()));
            }
            else
            {
                temp_triangle = Triangle_2(*iter, new_point, *(iter + 1));
                seg = Segment_2(*iter, *(iter + 1));
            }

            if(this->visible(seg, new_point))
            {
                double temp_area = std::abs(temp_triangle.area());
                if(area > temp_area)
                {
                    NewPolygon = Polygon_2(Real_Polygon);
                    area = temp_area;
                    NewPolygon.insert((iter - Real_Polygon.begin() + 1) + NewPolygon.begin(), new_point);
                }
            }
        }

        Real_Polygon = NewPolygon;
    }
    //if we try to find the max area polygon
    else if(how_to_remove_edge == '3')
    {
        double area = 0.0;
        Polygon_2 NewPolygon(Real_Polygon);

        for(typename Polygon_2::Vertices::iterator iter = lower_limit_iter ; iter <= uper_limit_iter ; iter++)
        {
            Triangle_2 temp_triangle;
            Segment_2 seg;

            if(iter == Real_Polygon.vertices_end() - 1)
            {
                temp_triangle = Triangle_2(*iter, new_point, *(Real_Polygon.vertices_begin()));
                seg = Segment_2(*iter, *(Real_Polygon.vertices_begin()));
            }
            else
            {
                temp_triangle = Triangle_2(*iter, new_point, *(iter + 1));
                seg = Segment_2(*iter, *(iter + 1));
            }

            if(this->visible(seg, new_point))
            {
                double temp_area = std::abs(temp_triangle.area());
                if(area < temp_area)
                {
                    NewPolygon = Polygon_2(Real_Polygon);
                    area = temp_area;
                    NewPolygon.insert((iter - Real_Polygon.begin() + 1) + NewPolygon.begin(), new_point);
                }
            }
        }

        Real_Polygon = NewPolygon;
    }
}

template<class Kernel> bool Incremental<Kernel>::visible(const Segment_2 seg, const Point_2 new_point)
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

template<class Kernel> bool Incremental<Kernel>::red_visible(const Segment_2 seg, const Point_2 new_point) 
{
    CGAL::Orientation point_orientation = CGAL::orientation(new_point, seg[0], seg[1]);//find the relative position between the new point and the line created by the edge
    CGAL::Orientation polygon_orientation;

    //if point on top of the line presented by segment, then the semgment is not visible by the point
    if (point_orientation == CGAL::COLLINEAR) 
    {
        return false;
    }

    for (Point_2 vertex : Convex_Hull_Polygon.vertices()) 
    {
        if (seg[0] == vertex || seg[1] == vertex)
        {
            continue;
        }
        //find orientation between the line and the rest of the polygon
        polygon_orientation = CGAL::orientation(vertex, seg[0], seg[1]);
        break;
    }
    //if point and polygon have different orientations to the semgent, then edge is visible by the point
    return (point_orientation != polygon_orientation);

}

#endif