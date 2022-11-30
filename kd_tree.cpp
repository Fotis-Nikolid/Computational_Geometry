#ifndef KDTREE_CPP
#define KDTREE_CPP
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "kd_tree.h"


template<class Kernel> kdTree<Kernel>::Node Node(Point_2 p, Node* left, Node* right)
{
    this->point = p;
    this->left = left;
    this->right = right;
}

template<class Kernel> kdTree<Kernel>::Node ~Node()
{
    delete left;
    delete right;
}

template<class Kernel> compare_x(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2)
{
    return p1.x() < p2.x();
}
template<class Kernel> compare_y(const CGAL::Point_2<Kernel> p1, const CGAL::Point_2<Kernel> p2)
{
    return p1.y() < p2.y();
}

template<class Kernel> Pvector::iterator kdTree<Kernel>::find_median(Pvector::iterator begin, Pvector::iterator end, int x_or_y)
{
    if(x_or_y == 0)
    {
        std::sort(begin, end, compare_x<Kernel>);
    }
    else if(x_or_y == 1)
    {
        std::sort(begin, end, compare_y<Kernel>);
    }

    return ((end - begin)/2) + begin;
}

template<class Kernel> Node* kdTree<Kernel>::insert(Pvector::iterator begin, Pvector::iterator end, int depth = 0)
{
    if(begin >= end)
    {
        return nullptr;
    }

    Pvector::iterator it = find_median(begin, end, depth % 2);

    return new Node(*it, this->insert(begin, it, depth + 1), this->insert(it + 1, end, depth + 1));
}

template<class Kernel> kdTree<Kernel>::kdTree(Pvector points)
{
    Pvector v(points);

    this->root = this->insert(v.begin(), v.end(), 0);
}

template<class Kernel> Pvector kdTree<Kernel>::find_points_inside_bounds(int upper_x, int upper_y, int lower_x , int lower_y)
{
    
}