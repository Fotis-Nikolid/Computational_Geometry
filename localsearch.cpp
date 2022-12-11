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
#include "hull.h"

template <class Kernel>
bool comp_min(const CGAL::Polygon_2<Kernel> &p1, const CGAL::Polygon_2<Kernel> &p2)
{
    // p1 old best p2 temp
    return std::abs(p2.area()) < std::abs(p1.area());
}

template <class Kernel>
bool comp_max(const CGAL::Polygon_2<Kernel> &p1, const CGAL::Polygon_2<Kernel> &p2)
{
    return std::abs(p2.area()) > std::abs(p1.area());
}

template <class Kernel>
LocalSearch<Kernel>::LocalSearch(const std::vector<Point_2> &Points, const std::string &min_or_max, const std::string &InitAlgorithm) : InitFailed(false)
{
    Hull<Kernel> hull; 
    if (min_or_max == "min")
    {
        this->compare = comp_min;
        if(InitAlgorithm == "incremental")
        {
            Polygon = Incremental<Kernel>(Points, "1a", '2').getPolygon();
        }
        else
        {
            if (hull.solve(Polygon, std::list<Point_2>(Points.begin(),Points.end()), '2')==0.0) {
                InitFailed=true;
            }
        }
    }
    else if (min_or_max == "max")
    {
        this->compare = comp_max;
        if(InitAlgorithm == "incremental")
        {
            Polygon = Incremental<Kernel>(Points, "1a", '3').getPolygon();
        }
        else
        {
            if (hull.solve(Polygon,std::list<Point_2>(Points.begin(),Points.end()), '3')==0.0) {
                InitFailed=true;
            }
        }
    }

    std::srand(std::time(nullptr));
}

template <class Kernel>
bool LocalSearch<Kernel>::InitializationFailed()
{
    return InitFailed;
}

template <class Kernel>
bool LocalSearch<Kernel>::MinimizePolygon(const int L, const double threshold, const int K)
{
    if(InitFailed) return false;

    return this->solve(L, threshold, K);
}

template <class Kernel>
bool LocalSearch<Kernel>::solve(const int L, const double threshold, const int K)
{
    if(L == 0)
        return false;
    
    bool solved = false;
    double diff = 0.0;
    Polygon_2 BestPolygon;

    if (K == 0 || K >= Polygon.edges().size() - 1)
    {
        Polygon_2 temp;
        do
        {
            diff = 0.0;
            for (int edge_index = 0 ; edge_index < Polygon.edges().size() ; edge_index++)
            {
                temp = Polygon;
                // swap the random edge with the best L
                double new_diff = ReplaceEdgeWithBest_L(&temp, edge_index, L);

                //if we found an L that has better area than the current Polygon
                if (new_diff > diff)
                {
                    BestPolygon = temp;
                    diff = new_diff;
                    solved = true;
                }
            }

            if (diff > 0.0)
                Polygon = BestPolygon;

        } while (diff > threshold);
    }
    else
    {
        Polygon_2 temp;
        do
        {
            diff = 0.0;
            std::unordered_map<int, bool>random_edge_indexes; //map that holds the random edge picks so we dont repick them
            //find K random edges and swap them with L
            for(int i = 0 ; i < K ; i++)
            {
                int pick = std::rand()%Polygon.edges().size();

                //if the segment has been picked again try another segment
                while(random_edge_indexes.find(pick) != random_edge_indexes.end())
                {
                    pick = std::rand()%Polygon.edges().size();
                }

                temp = Polygon;
                //swap the random edge with the best L
                double new_diff = ReplaceEdgeWithBest_L(&temp, pick, L);

                //if we found an L that has better area than the current Polygon
                if(new_diff > diff)
                {
                    solved = true;
                    BestPolygon = temp;
                    diff = new_diff;
                }
                random_edge_indexes[pick] = true;
            }

            if (diff > 0.0)
                Polygon = BestPolygon;

        }while(diff > threshold);
    }

    return solved;
}

//this function is helper function only used inside ReplaceEdgeWithBest_L
template <class Kernel>
bool LinesIntersect(const CGAL::Segment_2<Kernel>& seg1, const CGAL::Segment_2<Kernel>& seg2)
{
    auto intersect_p = intersection(seg1, seg2);
    if (intersect_p)
    {
        if (CGAL::Point_2<Kernel> *p = boost::get<CGAL::Point_2<Kernel>>(&*intersect_p))
        {
            if (*p == seg1[1] && seg2[1] == seg1[1])
                return false;
            else if(*p == seg1[0] && seg2[0] == seg1[0])
                return false;
            else
                return true;
        }
        else
        {
            return true;
        }
    }
    return false;
}

template <class Kernel>
double LocalSearch<Kernel>::ReplaceEdgeWithBest_L(Polygon_2 *BestPol, const int edge_destroy, const int L)
{
    // edge_destoy points to the location of the Segment_2[0] , where the segment is the edge we will destroy and put L
    const int dist = Polygon.vertices_end() - Polygon.vertices_begin();

    Polygon_2 temp;

    double diff = 0.0;

    for(int l = L ; l > 0 ; l--)
    {
        if (l >= dist - 1)
        {
            l = dist - 1;
            continue;
        }

        for (int Lstart = 0; Lstart < dist; Lstart++)
        {
            int before_Lstart = Lstart - 1;
            // if Lstart was the first polygon vertex the set before_Lstart to the last polygon vertex
            if (before_Lstart == -1)
            {
                before_Lstart = dist - 1;
            }

            int edge_end = edge_destroy + 1;
            // if we will destoy the last edge of the polygon then the Segment[1] = first polygon vertex
            if (edge_end == dist)
            {
                edge_end = 0;
            }

            int Lend = Lstart + (l - 1);
            // if L Lends after the last Polygon vertex
            if (Lend >= dist)
            {
                // set new Lend
                // if we asume remaining elements = total elements - (the elements from Lstart to the last vertex)
                // then the Lend will be in the first vertex + remaining elements
                Lend = Lend - dist;
                // if edge to destory is included in L then stop searching you have found all avaible L
                if (edge_destroy <= Lend)
                {
                    break;
                }
            }
            else
            {
                // if edge to destory is included in L then skip all the L that includes it
                if ((edge_destroy <= Lend && edge_destroy >= Lstart) || edge_end == Lstart)
                {
                    Lstart = edge_destroy + 1;
                    continue;
                }
            }

            int after_Lend = Lend + 1;
            // if Lend was the last polygon vertex then after_Lend to the first polygon vertex
            if (after_Lend == dist)
            {
                after_Lend = 0;
            }

            //3 new edges will be created , as we put L in the place of the edge we choose to destroy
            Segment_2 NewEdge1(*(Polygon.vertices_begin() + edge_destroy), *(Polygon.vertices_begin() + Lend));
            Segment_2 NewEdge2(*(Polygon.vertices_begin() + edge_end), *(Polygon.vertices_begin() + Lstart));
            Segment_2 NewEdge3(*(Polygon.vertices_begin() + before_Lstart), *(Polygon.vertices_begin() + after_Lend));
            //if those 3 edges intersect with one another then find onother L
            if(LinesIntersect(NewEdge1, NewEdge2)
            || LinesIntersect(NewEdge1, NewEdge3)
            || LinesIntersect(NewEdge2, NewEdge3)
            )
            {
                continue;
            }

            std::unordered_map<Segment_2, bool> map;
            map[(Segment_2)*(Polygon.edges_begin() + before_Lstart)] = true;
            map[(Segment_2)*(Polygon.edges_begin() + Lend)] = true;
            map[(Segment_2)*(Polygon.edges_begin() + edge_destroy)] = true;

            // check if the 3 new edges are visible
            if (!this->visible_points(NewEdge3[0], NewEdge3[1], map) 
                || !this->visible_points(NewEdge1[0], NewEdge1[1], map) 
                || !this->visible_points(NewEdge2[0], NewEdge2[1], map)
            )
            {
                continue;
            }

            temp = Polygon;
            // swap L potition
            RelocateEdges(temp, Lstart, Lend, edge_destroy, l);

            // check if new polygon is better
            if (this->compare(*BestPol, temp))
            {
                *BestPol = temp;
                diff = abs(abs(Polygon.area()) - abs(BestPol->area()));
            }
        }
    }

    return diff;
}


template <class Kernel>
void LocalSearch<Kernel>::RelocateEdges(Polygon_2 &Pol, int Lstart, int Lend, int edge_destroy, int L)
{
    // if start is after Lend then L Lends after the polygon start
    if (Lend < Lstart)
    {
        Pol.erase(Pol.vertices_begin() + Lstart, Pol.vertices_end());
        Pol.erase(Pol.vertices_begin(), Pol.vertices_begin() + Lend + 1);

        // after we remove L points then the edge destroy location moves Lend + 1 seats to the left
        edge_destroy -= Lend + 1;

        for(typename Polygon_2::Vertex_iterator it = Polygon.vertices_begin() + Lstart ; it < Polygon.vertices_end() ; it++)
        {
            Pol.insert(Pol.vertices_begin() + edge_destroy + 1, *it);
        }

        for(typename Polygon_2::Vertex_iterator it = Polygon.vertices_begin() ; it <= Polygon.vertices_begin() + Lend ; it++)
        {
            Pol.insert(Pol.vertices_begin() + edge_destroy + 1, *it);
        }
    }
    // else if L Lends behind polygon Lend
    else
    {
        Pol.erase(Pol.vertices_begin() + Lstart, Pol.vertices_begin() + Lend + 1);

        // if edge destroy is after the L then we need to move it
        if (edge_destroy > Lend)
        {
            // move L seats
            edge_destroy -= L;
        }

        for(typename Polygon_2::Vertex_iterator it = Polygon.vertices_begin() + Lstart ; it <= Polygon.vertices_begin() + Lend ; it++)
        {
            Pol.insert(Pol.vertices_begin() + edge_destroy + 1, *it);
        }

    }
}

template <class Kernel>
float LocalSearch<Kernel>::getPolygonArea()
{
    return abs(Polygon.area());
}

template <class Kernel>
CGAL::Polygon_2<Kernel> LocalSearch<Kernel>::getPolygon()
{
    return Polygon;
}

template <class Kernel>
bool LocalSearch<Kernel>::visible_points(const Point_2 &p1, const Point_2 &p2, const std::unordered_map<Segment_2, bool>& excluded_edges)
{
    Segment_2 seg = Segment_2(p1, p2);
    for (Segment_2 PolygonEdge : Polygon.edges())
    {

        if(excluded_edges.find(PolygonEdge) != excluded_edges.end())
            continue;

        auto intersect_p = intersection(PolygonEdge, seg);
        if (intersect_p)
        {
            if (Point_2 *p = boost::get<Point_2>(&*intersect_p))
            {
                if (*p != p1 && *p != p2)
                    return false;
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