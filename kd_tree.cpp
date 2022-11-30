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

template<class Kernel> Pvector::iterator kdTree<Kernel>::find_median(Pvector::iterator begin, Pvector::iterator end)
{
    
}