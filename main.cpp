#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include "polygon.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point_2;

int main(int argc, char *argv[])
{
  std::ifstream infile;
  std::ofstream outfile;
  std::string algorithm;
  std::string edge_selection;
  std::string sorting;

  for(int i = 1 ; i < argc ; i++)
  {
    if(argv[i][0] == '-' && (i + 1) < argc)
    {
      std::string arg = argv[i];
      std::string opt = argv[i + 1];

      if(arg == "-i")
      {
        infile.open(opt);
      }
      else if(arg == "-o")
      {
        outfile.open(opt);
      }
      else if(arg == "-algorithm")
      {
        algorithm=opt;
      }
      else if(arg == "-edge_selection")
      {
        edge_selection=opt;
      }
      else
      {
        std::cout << "Option " + arg + " is not avaible" << std::endl;
        return -1;
      }
    }
  }

  if (!infile.is_open())
  {
    std::cerr << "There was a problem opening the input file!\n";
    return -1;
  }

  if(!outfile.is_open())
  {
    std::cerr << "There was a problem opening the output file!\n";
    return -1;
  }

  std::list<Point_2> l_points;
  std::vector<Point_2> v_points;
  std::string line;
  while (getline(infile,line)) {
      double x,y;
      int row;
      std::stringstream sline(line);
      sline >> row;
      sline >> x;
      sline >> y;

      v_points.push_back(Point_2(x,y));
      l_points.push_back(Point_2(x,y));
  }

  Polygon<K> poly;
  if (algorithm=="hull") {
    poly.Hull_Based(v_points,edge_selection);
    poly.Size();
  }
  else {

  }

  
  return 0;
}
