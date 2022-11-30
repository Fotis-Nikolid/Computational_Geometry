#ifndef KDTREE_H
#define KDTREE_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <vector>

template<class Kernel> class kdTree
{
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef std::vector<Point_2> Pvector;

    private:
        class Node
        {
            public:
                Point_2 point;
                Node* left;
                Node* right;
                Node(Point_2, Node*, Node*);
                ~Node();
        };

        Node* root;

        Pvector::iterator find_median(Pvector::iterator, Pvector::iterator);
        Node insert(Pvector::iterator, Pvector::iterator, int);
    
    public:
        kdTree(Pvector);
        ~kdTree();
        Pvector find_points_inside_bounds(int, int, int, int);//upper x , y lower x , y (in this order)
};

#include "kd_tree.cpp"

#endif