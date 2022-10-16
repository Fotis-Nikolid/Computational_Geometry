/*#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>

template<class Kernel>
class Incremental {
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef Kernel::Point_2  Point_2;
    typedef Kernel::Segment_2 Segment;

    private: 
        //members
        Polygon_2 Convex_Hull;
        float Area;

        //internal helper functions 
        void Initialize(Polygon_2&);
        void Sort_Axis(std::vector<Point_2>,std::string);
        void Red_Edges(Polygon_2,Point_2);
        bool is_visible(Point_2,Point_2,Point_2);
    public:
        float solve(Polygon_2&,std::vector<Point_2>&,std::string);
};
*/