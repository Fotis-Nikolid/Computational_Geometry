/*#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::vector<Point_2> Points;

int main(int argc, char *argv[])
{
  std::ifstream infile;
  std::ofstream outfile;
  int algorithm {0};

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
        if(opt == "incremental")
        {
          algorithm = 1;
        }
        else if(opt == "convex_hull")
        {
          algorithm = 2;
        }
        else
        {
          std::cout << "Algorithm " + opt + " is not implemented" << std::endl;
          return -1;
        }
      }
      else if(arg == "-edge_selection")
      {
        //do something here
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

  Points pv;
  std::string line;

  if(!std::getline(infile,line))
  {
    std::cerr << "Something worng with the input second line" << std::endl;
    return -1;
  }

  // get polygon area
  if(std::getline(infile,line))
  {
    
  }
  else
  {
    std::cerr << "Something worng with the input second line" << std::endl;
    return -1;
  }

  //get the points
  while(std::getline(infile,line))
  {
    double x,y;
    std::stringstream sline(line);
    sline >> x;
    sline >> y;

    pv.push_back(Point_2(x,y));
  }

  for (std::vector<Point_2>::iterator it = pv.begin() ; it != pv.end(); ++it)
  {
    outfile << *it  << std::endl;
  }

  return 0;
}
*/