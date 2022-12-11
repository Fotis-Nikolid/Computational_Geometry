#if 1
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <chrono>
#include <regex>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include "polygon.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_2<Kernel> Point_2;
typedef CGAL::Segment_2<Kernel> Segment_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

int main(int argc, char *argv[])
{
  std::ifstream infile;
  std::ofstream outfile;
  std::string algorithm("simulated_annealing");
  std::string criteria("max");
  std::string initialization("convex_hull");
  int attempts=1;
  int L = 1;
  int K = 0;
  double threshold = 0.0;
  std::string step_choice("subdivision");

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
      else if(arg == "-initialization")
      {
        initialization=opt;
      }
      else if(arg == "-annealing")
      {
        step_choice=opt;
      }
      else if(arg == "-L")
      {
        L=atoi(opt.c_str());
      }
      else if(arg == "-attempts")
      {
        attempts=atoi(opt.c_str());
      }
      else if(arg == "-min")
      {
        criteria="min";
      }
      else if(arg == "-max")
      {
        criteria="max";
      }
      else if(arg == "-threshold")
      {
        threshold=atol(opt.c_str());
      }
      else if(arg == "-K")
      {
        K = atoi(opt.c_str());
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
  getline(infile,line);
  getline(infile,line);
  double convex_hull_area;
  std::smatch match_obj;
  if (std::regex_search(line,match_obj,std::regex("\"([0-9]*)\""))) {
      convex_hull_area=atof(match_obj[1].str().c_str());
  }
  while (getline(infile,line)) {
      double x,y;
      int row;
      std::stringstream sline(line);
      sline >> row;
      sline >> x;
      sline >> y;

      v_points.push_back(Point_2(x,y));
  }

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  Polygon<Kernel> poly(v_points,algorithm,criteria,step_choice,initialization,L,threshold,K,attempts);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  if (!poly.Success()) {
    std::cout<<"Polygon Creation Failed..."<<std::endl;
    return -1;
  }
  Polygon_2 real_poly(poly.get_Polygon());
  if (poly.Simple() && poly.Size() == v_points.size()) {
    std::cout<<"Correct"<<std::endl;
  }
  else if (poly.Size() != v_points.size()) {
    std::cout<<"Error: Wrong Number of points "<<poly.Size()<<std::endl;
    return -1;
  }
  else {
    std::cout<<"Error: Not simple"<<std::endl;
    return -1;
  }

  outfile<<"Polygonization"<<std::endl;
  for (const Point_2 p: real_poly.vertices()) {
    outfile<<p<<std::endl;
  }
  for (const Segment_2 edge: real_poly.edges()) {
    outfile<<edge<<std::endl;
  }
  outfile<<"Algorithm: "<<algorithm<<"_"<<criteria<<std::endl;

  outfile<<std::endl;

  double initial_area=poly.Init_Area();
  #if 1 //here we find the initial are for the case of subdivision, where it is not possible to extract initial area from the sum of initial areas of the sum polygons
    std::cout<<"Optimal Polygon Creation Finished, calculating Initial Area....(This may take a while)"<<std::endl;
    if (step_choice=="subdivision") {
      char crit=(criteria=="max")?('3'):('2');
      if (initialization=="incremental") {
        Incremental<Kernel> inc(v_points,"1a",crit);
        initial_area=inc.getPolygonArea();
      }
      else if (initialization=="convex_hull") {
        Hull<Kernel> hull;
        std::list<Point_2> points(v_points.begin(),v_points.end());
        Polygon_2 t_poly;
        initial_area=hull.solve(t_poly,points,crit,NULL,NULL);
      }
    }
  #endif
  outfile<<"area: "<<poly.Area()<<std::endl;
  outfile<<"area_initial: "<<initial_area<<std::endl;
  outfile<<"ratio: "<<poly.Area()/convex_hull_area<<std::endl;
  outfile<<"ratio_initial: "<<initial_area/convex_hull_area<<std::endl;
  outfile<<"construction time: "<<std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<<std::endl;

  return 0;
}

#endif