#ifndef LOCALSEARCH_CPP
#define LOCALSEARCH_CPP
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "localsearch.h"
#include "incremental.h"


template<class Kernel> bool comp_min(const CGAL::Polygon_2<Kernel>& p1, const CGAL::Polygon_2<Kernel>& p2)
{
    //p1 old best p2 temp
    return std::abs(p2.area()) < std::abs(p1.area());
}

template<class Kernel> bool comp_max(const CGAL::Polygon_2<Kernel>& p1, const CGAL::Polygon_2<Kernel>& p2)
{
    return std::abs(p2.area()) > std::abs(p1.area());
}

template<class Kernel> LocalSearch<Kernel>::LocalSearch(const std::vector<Point_2>& Points, const std::string& min_or_max, const int L, const int K, const double threshold) : Polygon(Incremental<Kernel>(Points, "1a", '2').getPolygon())
{
    if(min_or_max == "min")
    {
        this->compare = comp_min;
    }
    else if(min_or_max == "max")
    {
        this->compare = comp_max;
    }

    std::cout << getPolygonArea() << std::endl;

    std::srand(std::time(nullptr));

    for(int k = K ; k > 0 ; k--)
    {
        bool solved = this->solve_specific_K(L, k, threshold);
        //if a better polygon can not be found with this L try L-1 until one is found or L goes to 0
        for(int l = L - 1 ; l > 0 && solved == false ; l--)
        {
            solved = this->solve_specific_K(l, k, threshold);
        }
    }
}

template<class Kernel> bool LocalSearch<Kernel>::solve_specific_K(const int L, const int K , const double threshold)
{
    bool solved = false;
    double diff;

    do
    {
        diff = 0.0;
        int old_choice = -1;
        for(int i = 0 ; i < K ; i++)
        {
            //pick random edge
            if(temp.size() == 0)
            {
                break;
            }
            
            std::swap(temp[vertices_pot.size() - 1], temp[std::rand()%temp.size()]);

            double sr = swap_L_with_edge(temp[temp.size() - 1], L);

            if(sr != 0.0)
            {
                solved = true;
                diff = sr;
            }

            vertices_pot.pop_back();
        }
    }while(diff >= threshold);

    return solved;
}

template<class Kernel> void LocalSearch<Kernel>::relocate_edges(Polygon_2& Pol, int Lstart, int Lend, int edge_destroy, int L)
{

    //if start is after end then L ends after the polygon start
    if(Lend < Lstart)
    {
        //after we remove L points then the edge destroy location moves end + 1 seats to the left
        edge_destroy -= Lend + 1;
    }
    //else if L ends behind polygon end
    else
    {
        //if edge destroy is after the L then we need to move it
        if(edge_destroy > Lend)
        {
            //move L seats
            edge_destroy -= L;
        } 
    }

    Pol.erase(Pol.vertices_begin() + Lstart, Pol.vertices_begin() + Lend + 1);

    Pol.insert(Pol.vertices_begin() + edge_destroy + 1, Polygon.vertices_begin() + Lstart, Polygon.vertices_begin() + Lend + 1);

}

template<class Kernel> double LocalSearch<Kernel>::swap_L_with_edge(const int edge_destroy, const int L)
{
    const int dist = Polygon.vertices_end() - Polygon.vertices_begin();

    if(L >= dist)
    {
        return 0.0;
    }

    double diff = 0.0;

    Polygon_2 BestPol(Polygon);

    for(int i = 0 ; i < dist ; i++)
    {
        int before_i  = i - 1;

        if(before_i == -1)
        {
            before_i = dist - 1;
        }

        int end = i + (L - 1);
        if(end >= dist)
        {
            end = end - dist;
            //debuging
            if(end < 0)
            {
                std::cerr << "End less that 0 in swap" << std::endl;
            }
            //if edge to destory is included in L then stop searching you have found all avaible L 
            if(edge_destroy >= i || edge_destroy <= end)
            {
                break;
            }
        }
        else
        {
            //if edge to destory is included in L then skip all the L that includes it 
            if(edge_destroy <= end && edge_destroy >= i)
            {
                i = edge_destroy + 1;
                continue;
            }
        }

        int after_end = end + 1;
        if(after_end == dist)
        {
            after_end = 0;
        }

        int edge_end = edge_destroy + 1;
        if(edge_end == dist)
        {
            edge_end = 0;
        }

        //check if 3 new edges are visible
        if(this->visible_points(*(Polygon.vertices_begin() + before_i), *(Polygon.vertices_begin() + after_end)) 
           && this->visible_points(*(Polygon.vertices_begin() + i), *(Polygon.vertices_begin() + edge_destroy)) 
           && this->visible_points(*(Polygon.vertices_begin() + end), *(Polygon.vertices_begin() + edge_end))
          )
        {
            Polygon_2 temp(Polygon);
            //swap L potition
            relocate_edges(temp, i, end, edge_destroy, L);
            //check if new polygon is better
            if(this->compare(BestPol, temp))
            { 
                BestPol = temp;
                diff = abs(abs(Polygon.area()) - abs(BestPol.area()));
            }
        }
    }

    if(diff != 0.0)
        Polygon = BestPol;
    
    return diff;
}

template<class Kernel> float LocalSearch<Kernel>::getPolygonArea()
{
    return abs(Polygon.area());
}

template<class Kernel> CGAL::Polygon_2<Kernel> LocalSearch<Kernel>::getPolygon()
{
    return Polygon;
}

template<class Kernel> bool LocalSearch<Kernel>::visible_points(const Point_2& p1, const Point_2& p2)
{
    Segment_2 seg = Segment_2(p1, p2);
    for(Segment_2 PolygonEdge : Polygon.edges())
    {
        
        auto intersect_p = intersection(PolygonEdge, seg);
        if(intersect_p)
        {
            if(Point_2* p=boost::get<Point_2 >(&*intersect_p))
            {
                if(*p != p1 && *p != p2) return false;
            }
            else
            {
                return false;
            }
        }
    }

    return true;
}

#endif