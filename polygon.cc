#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include "hull.h"
#include "polygon.h"


template<class Kernel>
void Polygon<Kernel>::Hull_Based(std::vector<Point_2> Points,std::string Edge_Selection) {
    Hull<Kernel> alg;
    std::list<Point_2> list(Points.begin(),Points.end());
    dt_Area=alg.solve(Polygon,list,*(Edge_Selection.c_str()));
}
template<class Kernel>
void Polygon<Kernel>::Increment_Based(std::vector<Point_2> Points,std::string Edge_Selection,std::string Sorting) {


}
template<class Kernel>
int Polygon<Kernel>::Size() {
    return Polygon.vertices().size();
}
template<class Kernel>
std::vector<CGAL::Point_2<Kernel>> Polygon<Kernel>::Points() {
    return Polygon.vertices();
}
template<class Kernel>
CGAL::Polygon_2<Kernel> Polygon<Kernel>::get_Polygon() {
    return this->Polygon;
}
template<class Kernel>
float Polygon<Kernel>::Area() {
    return dt_Area;
}
template<class Kernel>
float Polygon<Kernel>::Ratio() {
    
}
template<class Kernel>
bool Polygon<Kernel>::Simple() {
    return Polygon.is_simple();    
}



