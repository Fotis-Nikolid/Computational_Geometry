#if 1
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <list>
#include <regex>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include "localsearch.h"
#include "simulated_annealing.h"
#include <chrono>
#include <map>
#include <vector>

template<class Kernel> class Report {
  typedef CGAL::Point_2<Kernel> Point_2;
  typedef CGAL::Polygon_2<Kernel> Polygon_2;
  
  private:
    std::vector<std::string> Min_Max={"min","max"};
    std::vector<std::string> algorithms={"Hull+Subdv+Anneal(loc)","Hull+Subd+LocSearch","Hull+Anneal(loc)+LocSearch","Inc+Anneal(glob)+LocSearch"};
    std::map<int,std::vector<int>> min_score;
    std::map<int,std::vector<int>> max_score;
    std::map<int,std::vector<int>> min_bound;
    std::map<int,std::vector<int>> max_bound;

    void proccess_result(std::string criteria,int n_points,int alg_n,double ratio) {
      if (criteria=="min") {
        if (min_score.find(n_points)==min_score.end()) {
          min_score.find(n_points)->second.push_back(ratio);
        }
        else {
          min_score.find(n_points)->second.at(alg_n)+=ratio;
        }
        if (min_bound.find(n_points)==min_bound.end()) {
          min_bound.find(n_points)->second.push_back(ratio);
        }
        else if (ratio<min_bound.find(n_points)->second.at(alg_n)) {
          min_bound.find(n_points)->second.at(alg_n)=ratio;
        }
      }
      else {
        if (max_score.find(n_points)==max_score.end()) {
          max_score.find(n_points)->second.push_back(ratio);
        }
        else {
          max_score.find(n_points)->second.at(alg_n)+=ratio;
        }
        if (max_bound.find(n_points)==max_bound.end()) {
          max_bound.find(n_points)->second.push_back(ratio);
        }
        else if (ratio<max_bound.find(n_points)->second.at(alg_n)) {
          max_bound.find(n_points)->second.at(alg_n)=ratio;
        }
      }
    }
  public:
    void create_report(std::string dir_path,std::string outfile_path,bool preprocess) {
      DIR *dir;
      struct dirent *ent;
      std::ofstream outfile;
      outfile.open(outfile_path);

      //default parameters in case preproccessing is not performed
      int attempts=1;
      int annealing_L = 1000;
      int search_L=3;
      int K = 10;
      double threshold = 10.0;
      double initial_area;

      std::string file_name;
      if ((dir = opendir(dir_path.c_str())) != NULL) {//work on all files inside the directory(in random order) 
        while ((ent = readdir (dir)) != NULL) {
          file_name=ent->d_name;
          if (file_name=="." || file_name=="..") {
            continue;
          }
          
          std::string file_path=dir_path+"/"+file_name;
          std::ifstream file;
          file.open(file_path);
          if (!file.is_open()) {
            std::cerr << "There was a problem opening the input file!\n";
            return ;
          }
          std::string line;
          
          getline(file,line);
          getline(file,line);
          double convex_area;
          std::smatch match_obj;
          if (std::regex_search(line,match_obj,std::regex("\"([0-9]*)\""))) {//get area of convex hull
            convex_area=atof(match_obj[1].str().c_str());
          }
      
          std::vector<Point_2> v_points;
          while (getline(file,line)) {
            double x,y;
            int row;
            std::stringstream sline(line);
            sline >> row;
            sline >> x;
            sline >> y;

            v_points.push_back(Point_2(x,y));
          }
          long int n_points=v_points.size();//get number of points


          if (preprocess) {
            //do preproccessing and update optimal hyper-parameters 
          }        
          std::chrono::milliseconds cut_off;
          bool success;          
          //Hull+Subdivision+SimulatedAnnealing(local) 
          int alg_n=0;
          double area;
          Polygon_2 polygon_alt1,polygon_alt2,polygon_alt3,polygon_alt4;
          for (auto criteria:Min_Max) {
            Simulated_Annealing<Kernel> algorithm;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if ((success=algorithm.sub_division(polygon_alt1,v_points,criteria,"convex_hull",annealing_L,initial_area,cut_off))) { 
              area=abs(polygon_alt1.area());
            
              auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(start -std::chrono::system_clock::now());
              if (milliseconds<cut_off) {
                cut_off-=milliseconds;
                if ((success=algorithm.solve(polygon_alt1,v_points,criteria,"local","",annealing_L,1,initial_area,cut_off))) 
                  area=abs(polygon_alt1.area());
              }
            }
            double ratio;
            if (!success) {
              ratio=((criteria=="min")?1:0);
            }
            else {
              ratio=area/convex_area;
            }
            proccess_result(criteria,n_points,alg_n,ratio);
          }
          
          
          
          //Hull+Subdivision+LocalSearch
          alg_n++;
          for (auto criteria:Min_Max) {
            Simulated_Annealing<Kernel> algorithm;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if ((success=algorithm.sub_division(polygon_alt2,v_points,criteria,"convex_hull",annealing_L,initial_area,cut_off))) { 
              area=abs(polygon_alt2.area());
              auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(start -std::chrono::system_clock::now());
              if (milliseconds<cut_off) {
                cut_off-=milliseconds;
                if ((success=algorithm.solve(polygon_alt2,v_points,criteria,"local","convex_hull",annealing_L,1,initial_area,cut_off))) {
                  area=abs(polygon_alt2.area());

            
                  milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(start -std::chrono::system_clock::now());
                  if (milliseconds<cut_off) {
                    cut_off-=milliseconds;
                    LocalSearch<Kernel> loc(polygon_alt1,criteria);
                    if (!loc.MinimizePolygon(search_L, threshold,cut_off,K)) 
                      area = loc.getPolygonArea();
                  }
                }
              }
            }
            double ratio;
            if (!success) {
              ratio=(criteria=="min")?1:0;
            }
            else {
              ratio=area/convex_area;
            }
            proccess_result(criteria,n_points,alg_n,ratio);
          } 
          
          
          //Incremental+SimulatedAnnealing(global)+LocalSearch
          alg_n++;
          for (auto criteria:Min_Max) {
            Simulated_Annealing<Kernel> algorithm;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if (algorithm.solve(polygon_alt3,v_points,criteria,"global","incremental",annealing_L,1,initial_area,cut_off)) 
              area=abs(polygon_alt3.area());
            
            auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(start -std::chrono::system_clock::now());
            if (milliseconds<cut_off) {
              cut_off-=milliseconds;
              LocalSearch<Kernel> loc(polygon_alt3,criteria);
              if (!loc.MinimizePolygon(search_L, threshold,cut_off,K)) 
                area = loc.getPolygonArea();
            }
            double ratio;
            if (!success) {
              ratio=(criteria=="min")?1:0;
            }
            else {
              ratio=area/convex_area;
            }
            proccess_result(criteria,n_points,alg_n,ratio);
          }

          
          
          //Hull+SimulatedAnnealing(local)+LocalSearch
          alg_n++;
          for (auto criteria:Min_Max) {
            Simulated_Annealing<Kernel> algorithm;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if (algorithm.solve(polygon_alt4,v_points,criteria,"local","convex_hull",annealing_L,1,initial_area,cut_off)) 
              area=abs(polygon_alt4.area());
            auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(start-std::chrono::system_clock::now());
            if (milliseconds<cut_off) {
              cut_off-=milliseconds;
              LocalSearch<Kernel> loc(polygon_alt4,criteria);
              if (!loc.MinimizePolygon(search_L,threshold,cut_off,K)) 
                area = loc.getPolygonArea();
            }
            double ratio;
            if (!success) {
              ratio=(criteria=="min")?1:0;
            }
            else {
              ratio=area/convex_area;
            }
            proccess_result(criteria,n_points,alg_n,ratio);
          }
        }
        closedir (dir);
      } 
      else {        // could not open directory 
        perror ("");
        return;
      }
    }
};


int main(int argc,char* argv[]) {//receives starting file number, ending file number and a directory, and performs selected algorithm with all criteria for all files
  std::string dir_path;
  std::string outfile_path;
  bool preprocess=false;
  
  for(int i = 1 ; i < argc ; i++) {
    if(argv[i][0] == '-' && (i + 1) < argc) {
      std::string arg = argv[i];
      std::string opt = argv[i + 1];

      if (arg == "-i") {
        dir_path=opt;
      }
      else if (arg == "-o") {
        outfile_path=opt;
      }
      else if (arg=="-preprocess") {
        preprocess=true;
      }
      else {
        std::cout << "Option " + arg + " is not avaible" << std::endl;
        return -1;
      }
    }
  }
  Report<CGAL::Exact_predicates_inexact_constructions_kernel> rp;
  rp.create_report(dir_path,outfile_path,preprocess);
  return 0;
}
#endif