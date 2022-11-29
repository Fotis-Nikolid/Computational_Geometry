#ifndef LOCALSEARCH_CPP
#define LOCALSEARCH_CPP
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "LocalSearch.h"
#include "incremental.h"


template<class Kernel> bool comp_min(const CGAL::Polygon_2<Kernel> p1, const CGAL::Polygon_2<Kernel> p2)
{
    return std::abs(p2.area()) < std::abs(p1.area());
}

template<class Kernel> bool comp_max(const CGAL::Polygon_2<Kernel> p1, const CGAL::Polygon_2<Kernel> p2)
{
    return std::abs(p2.area()) > std::abs(p1.area());
}

template<class Kernel> LocalSearch<Kernel>::LocalSearch(const std::vector<Point_2> Points, const std::string min_or_max, const int L, const int K, const double threshold)
{
    if(min_or_max == "min")
    {
        this->compare = comp_min;
    }
    else if(min_or_max == "max")
    {
        this->compare = comp_max;
    }

    this->Polygon = Polygon_2(Incremental(Points, "1a", '1').getPolygon());

    std::srand(std::time(nullptr));

    for(int k = K ; k > 0 ; k--)
    {
        bool solved = solve_specific_K(min_or_max, L, k, threshold);
        //if a better polygon can not be found with this L try L-1 until one is found or L goes to 0
        for(int l = L - 1 ; l > 0 && solved == false ; l--)
        {
            solved = solve_specific_K(min_or_max, l, k, threshold);
        }
    }
}

template<class Kernel> bool LocalSearch<Kernel>::solve_specific_K(const int L, const int K , const double threshold)
{
    bool solved = false;
    double diff = 0.0;

    std::vector<int> vertices_pot(Polygon.vertices().size(), 0);

    for(int i = 0 ; i < Polygon.vertices().size() ; i++)
    {
        vertices_pot[i] = i;
    }

    do
    {
        for(int i = 0 ; i < K ; i++)
        {
            //pick random edge
            std::swap(vertices_pot[vertices_pot.size() - 1], vertices_pot[std::rand()%vertices_pot.size()]);

            double sr = swap_L_with_edge(vertices_pot[vertices_pot.size() - 1], L);

            if(sr != 0.0)
            {
                solved = true;
                diff = sr;
            }

            vertices_pot.pop_back();
        }
    }while(diff >= threshold)

    return solved;
}

template<class Kernel> void relocate_edges(CGAL::Polygon_2<Kernel>& Pol, int start, int end, int edge_destroy)
{
    return 0.0;
}

template<class Kernel> double LocalSearch<Kernel>::swap_L_with_edge(const int edge_destroy, const int L)
{
    Polygon_2 temp(Polygon);
    Polygon_2 BestPol(Polygon);
    double diff = 0.0;
    const int dist = Polygon.vertices().end() - Polygon.vertices().begin();

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

        if(this->visible_points(*(Polygon.vertices().begin() + before_i), *(Polygon.vertices().begin() + after_end)) && this->visible_points(*(Polygon.vertices().begin() + i), *(Polygon.vertices().begin() + edge_destroy)) && this->visible_points(*(Polygon.vertices().begin() + after_end), *(Polygon.vertices().begin() + edge_end)))
        {
            relocate_edges(temp, i, end, edge_destroy);

            if(this->compare(BestPol, temp))
            {
                BestPol = Polygon_2(temp);
                diff = abs(abs(Polygon.area()) - abs(BestPol.area()));
            }
        }
    }

    if(diff != 0.0)
        Polygon = BestPol;
    
    return diff;
}