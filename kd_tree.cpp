#ifndef KDTREE_CPP
#define KDTREE_CPP
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "kd_tree.h"

template<class Kernel> class kdTree<Kernel>::Node
{
    public:
        Point_2 point;
        Node* left;
        Node* right;

        Node(const CGAL::Point_2<Kernel>& p, typename kdTree<Kernel>::Node* l, typename kdTree<Kernel>::Node* r) : point(p), left(l), right(r)
        {
        }
        ~Node()
        {
            delete left;
            delete right;
        }
};

template<class Kernel> bool compare_x(const CGAL::Point_2<Kernel>& p1, const CGAL::Point_2<Kernel>& p2)
{
    return p1.x() < p2.x();
}
template<class Kernel> bool compare_y(const CGAL::Point_2<Kernel>& p1, const CGAL::Point_2<Kernel>& p2)
{
    return p1.y() < p2.y();
}

template<class Kernel> typename std::vector<CGAL::Point_2<Kernel>>::iterator kdTree<Kernel>::find_median(PointsIterator begin, PointsIterator end, const int x_or_y)
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

template<class Kernel> typename kdTree<Kernel>::Node* kdTree<Kernel>::insert(PointsIterator begin, PointsIterator end, int depth)
{
    if(begin >= end)
    {
        return nullptr;
    }

    PointsIterator it = find_median(begin, end, depth % 2);

    return new Node(*it, this->insert(begin, it, depth + 1), this->insert(it + 1, end, depth + 1));
}

template<class Kernel> kdTree<Kernel>::kdTree(PointsVector points)
{
    this->root = this->insert(points.begin(), points.end(), 0);
}

template<class Kernel> void kdTree<Kernel>::points_inside_bounds(PointsVector& points, Node *traverse, const int upper_x, const int upper_y, const int lower_x , const int lower_y, int depth)
{
    if(traverse == nullptr)
        return;
    
    int x = traverse->point.x();
    int y = traverse->point.y();

    if(depth % 2 == 0)
    {
        if(x < lower_x)
        {
            this->points_inside_bounds(points, traverse->right, upper_x, upper_y, lower_x, lower_y, depth + 1);
        }
        else if(x > upper_x)
        {
            this->points_inside_bounds(points, traverse->left, upper_x, upper_y, lower_x, lower_y, depth + 1);
        }
        else
        {
            if(y <= upper_y && y >= lower_y)
            {
                points.push_back(traverse->point);
            }
            this->points_inside_bounds(points, traverse->left, upper_x, upper_y, lower_x, lower_y, depth + 1);
            this->points_inside_bounds(points, traverse->right, upper_x, upper_y, lower_x, lower_y, depth + 1);
        }
    }
    else if(depth % 2 == 1)
    {
        if(y < lower_y)
        {
            this->points_inside_bounds(points, traverse->right, upper_x, upper_y, lower_x, lower_y, depth + 1);
        }
        else if(y > upper_y)
        {
            this->points_inside_bounds(points, traverse->left, upper_x, upper_y, lower_x, lower_y, depth + 1);
        }
        else
        {
            if(x <= upper_x && x >= lower_x)
            {
                points.push_back(traverse->point);
            }
            this->points_inside_bounds(points, traverse->left, upper_x, upper_y, lower_x, lower_y, depth + 1);
            this->points_inside_bounds(points, traverse->right, upper_x, upper_y, lower_x, lower_y, depth + 1);
        }
    }
}

template<class Kernel> void kdTree<Kernel>::find_points_inside_bounds(PointsVector& points, const int upper_x, const int upper_y, const int lower_x , const int lower_y)
{
    this->points_inside_bounds(points, root, upper_x, upper_y, lower_x, lower_y);
}

template<class Kernel> kdTree<Kernel>::~kdTree()
{
    delete root;
}

#endif