#ifndef LOCALSEARCH_CPP
#define LOCALSEARCH_CPP
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
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

template<class Kernel> LocalSearch<Kernel>::LocalSearch(const std::vector<Point_2>& Points, const std::string& min_or_max, const int L, const int K, const double threshold)
{
    if(min_or_max == "min")
    {
        this->compare = comp_min;
        Polygon = Incremental<Kernel>(Points, "1a", '2').getPolygon();
    }
    else if(min_or_max == "max")
    {
        this->compare = comp_max;
        Polygon = Incremental<Kernel>(Points, "1a", '3').getPolygon();
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

template<class Kernel> bool LocalSearch<Kernel>::solve_specific_K(const int L, const double threshold, const int K)
{
    bool solved = false;
    double diff = 0.0;
    Polygon_2 BestPolygon;

    if(K == 0 || K >= Polygon.edges().size() - 1)
    {
        Polygon_2 temp;
        do
        {
            diff = 0.0;
            for(Segment_2 edge : Polygon.edges())
            {
                temp = Polygon;
                //swap the random edge with the best L
                double sr = swap_L_with_edge(&temp, edge[0], L);

                if(sr > diff)
                {
                    BestPolygon = temp;
                    diff = sr;
                    solved = true;
                }
            }

            if(diff > 0.0)
                Polygon = BestPolygon

        }while(diff > threshold);
    }
    else
    {
        /*do
        {
            Polygon_2 BestPolygon(Polygon);
            diff = 0.0;
            std::unordered_map<Segment_2, bool>seg_map; //map that holds the random edge picks so we dont repick them
            //find K random edges and swap them with L
            for(int i = 0 ; i < K ; i++)
            {
                int pick = std::rand()%Polygon.edges().size();
                Segment_2 seg = (Segment_2)*(Polygon.edges_begin() + pick);
                
                //if the segment has been picked again try another segment
                while(seg_map.find(seg) != seg_map.end())
                {
                    pick = std::rand()%Polygon.edges().size();
                    seg = (Segment_2)*(Polygon.edges_begin() + pick);
                }

                //swap the random edge with the best L
                double sr = swap_L_with_edge(pick, L);

                //if we found an L that has better than the current Polygon
                if(sr != 0.0)
                {
                    solved = true;
                    diff = sr;
                }
                else
                {
                    seg_map[seg] = true;
                }
            }

        }while(diff >= threshold);*/
    }

    return solved;
}

template<class Kernel> void LocalSearch<Kernel>::relocate_edges(Polygon_2& Pol, int Lstart, int Lend, int edge_destroy, int L)
{
    //if start is after end then L ends after the polygon start
    if(Lend < Lstart)
    {
        Pol.erase(Pol.vertices_begin() + Lstart, Pol.vertices_end());
        Pol.erase(Pol.vertices_begin(), Pol.vertices_begin() + Lend + 1);
        
        //after we remove L points then the edge destroy location moves end + 1 seats to the left
        edge_destroy -= Lend + 1;

        Pol.insert(Pol.vertices_begin() + edge_destroy + 1, Polygon.vertices_begin() + Lstart, Polygon.vertices_end());

        //calculate the new potion to insert after the first insert
        int new_pot = edge_destroy + 1 + ((Polygon.vertices_end() - Polygon.vertices_begin()) - Lstart);

        Pol.insert(Pol.vertices_begin() + new_pot, Polygon.vertices_begin(), Polygon.vertices_begin() + Lend + 1);
    }
    //else if L ends behind polygon end
    else
    {
        Pol.erase(Pol.vertices_begin() + Lstart, Pol.vertices_begin() + Lend + 1);
        //if edge destroy is after the L then we need to move it
        if(edge_destroy > Lend)
        {
            //move L seats
            edge_destroy -= L;
        }

        Pol.insert(Pol.vertices_begin() + edge_destroy + 1, Polygon.vertices_begin() + Lstart, Polygon.vertices_begin() + Lend + 1); 
    }
}

template<class Kernel> double LocalSearch<Kernel>::swap_L_with_edge(Polygon_2* BestPol, const int edge_destroy, const int L)
{
    //edge_destoy points to the location of the Segment_2[0] , where the segment is the edge we will destroy and put L
    const int dist = Polygon.vertices_end() - Polygon.vertices_begin();

    if(L >= dist - 1)
    {
        return 0.0;
    }

    Polygon_2 temp;

    double diff = 0.0;

    for(int i = 0 ; i < dist ; i++)
    {
        int before_i  = i - 1;
        //if i was the first polygon vertex the set before_i to the last polygon vertex
        if(before_i == -1)
        {
            before_i = dist - 1;
        }

        int end = i + (L - 1);
        //if L ends after the last Polygon vertex
        if(end >= dist)
        {
            //set new end 
            //if we asume remaining elements = total elements - (the elements from i to the last vertex)
            //then the end will be in the first vertex + remaining elements
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
        //if end was the last polygon vertex then after_end to the first polygon vertex
        if(after_end == dist)
        {
            after_end = 0;
        }

        int edge_end = edge_destroy + 1;
        //if we will destoy the last edge of the polygon then the Segment[1] = first polygon vertex
        if(edge_end == dist)
        {
            edge_end = 0;
        }

        //check if the 3 new edges are visible
        if(this->visible_points(*(Polygon.vertices_begin() + before_i), *(Polygon.vertices_begin() + after_end)) 
           && this->visible_points(*(Polygon.vertices_begin() + i), *(Polygon.vertices_begin() + edge_destroy)) 
           && this->visible_points(*(Polygon.vertices_begin() + end), *(Polygon.vertices_begin() + edge_end))
          )
        {
            temp = Polygon;
            //swap L potition
            relocate_edges(temp, i, end, edge_destroy, L);
            
            //check if new polygon is better
            if(this->compare(*BestPol, temp))
            {
                *BestPol = temp;
                diff = abs(abs(Polygon.area()) - abs(BestPol->area()));
            }
        }
    }
    
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