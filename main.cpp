#if 1
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <dirent.h>
#include <list>
#include <regex>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include "localsearch.h"
#include "simulated_annealing.h"
#include "hull.h"
#include "incremental.h"
#include <chrono>
#include <map>
#include <vector>


template<class Kernel> class Report {
  typedef CGAL::Point_2<Kernel> Point_2;
  typedef CGAL::Polygon_2<Kernel> Polygon_2;
  
  private:
    std::vector<std::string> Min_Max={"min","max"};
    std::vector<std::string> algorithms={"Hull+Subd+Anneal(loc)","Hull+Subd+LocSearch","Hull+Anneal(loc)+LocSearch","Inc+Anneal(glob)+LocSearch"};
    std::map<int,std::vector<double>> min_score;
    std::map<int,std::vector<double>> max_score;
    std::map<int,std::vector<double>> min_bound;
    std::map<int,std::vector<double>> max_bound;

    void proccess_result(std::string criteria,int n_points,int alg_n,double ratio) {
      if (criteria=="min") {
        if (min_score.find(n_points)==min_score.end()) {
          min_score[n_points]={0,0,0,0};
          min_score.find(n_points)->second.at(alg_n)=ratio;
        }
        else {
          min_score.find(n_points)->second.at(alg_n)+=ratio;
        }
        if (min_bound.find(n_points)==min_bound.end()) {
          min_bound[n_points]={0,0,0,0};
          min_bound.find(n_points)->second.at(alg_n)=ratio;
        }
        else if (ratio>min_bound.find(n_points)->second.at(alg_n)) {
          min_bound.find(n_points)->second.at(alg_n)=ratio;
        }
      }
      else {
        if (max_score.find(n_points)==max_score.end()) {
          max_score[n_points]={0,0,0,0};
          max_score.find(n_points)->second.at(alg_n)=ratio;
        }
        else {
          max_score.find(n_points)->second.at(alg_n)+=ratio;
        }
        if (max_bound.find(n_points)==max_bound.end()) {
          max_bound[n_points]={1,1,1,1};
          max_bound.find(n_points)->second.at(alg_n)=ratio;
        }
        else if (ratio<max_bound.find(n_points)->second.at(alg_n)) {
          max_bound.find(n_points)->second.at(alg_n)=ratio;
        }
      }
    }
    void write_results(std::string outfile_path) {
      std::ofstream report;
      report.open(outfile_path);

      if (!report.is_open()) {
        std::cerr<<"Error in opening file for writing"<<std::endl;
        return ;
      }

      report << std::fixed << std::setprecision(6) << std::endl;
      report<<"      ||";
      for (auto algorithm:algorithms) {
        int size=algorithm.length();
        int max=50;
        for (int i=1;i<std::floor((max-size)/2.0);i++) {
          report<<" ";
        }
        report<<algorithm;
        for (int i=1;i<std::round((max-size)/2.0);i++) {
          report<<" ";
        }
        report<<"||";
      }
      report<<std::endl;
      report<<"Size  ||";
      for (auto algorithm:algorithms) {
        report<<" min score | max score | min bound | max bound  ||";
      }
      report<<std::endl;
      for (auto set_scores:min_score) {
        int n_points=set_scores.first;
        report<<n_points;
        for (int i=0;i<(6-std::to_string(n_points).length());i++) {
          report<<" ";
        }
        report<<"||";
        for (int i=0;i<algorithms.size();i++) {
          double min_sc=min_score.find(n_points)->second.at(i);
          report<<" ";
          report<<min_sc;
          for (int k=0;k<(10-std::to_string(min_sc).length());k++) {
            report<<" ";
          }
          report<<"|";
          double max_sc=max_score.find(n_points)->second.at(i);
          report<<" ";
          report<<max_sc;
          for (int k=0;k<(10-std::to_string(max_sc).length());k++) {
            report<<" ";
          }
          report<<"|";
          double min_bd=min_bound.find(n_points)->second.at(i);
          report<<" ";
          report<<min_bd;
          for (int k=0;k<(10-std::to_string(min_bd).length());k++) {
            report<<" ";
          }
          report<<"|";
          double max_bd=max_bound.find(n_points)->second.at(i);
          report<<" ";
          report<<max_bd;
          for (int k=0;k<(10-std::to_string(max_bd).length());k++) {
            report<<" ";
          }
          report<<" ||";
        }
        report<<std::endl;
      }
      report.close();
    }
  public:
    void evaluate(std::string dir_path,std::string outfile_path,bool preprocess) {
      DIR *dir;
      struct dirent *ent;
      //default parameters in case preprocessing is not performed
      int attempts=1;
      long int loc_L  =1000000;
      long int glob_L =100000;
      long int subd_L =100000;
      int search_L=3;
      int K=0;
      double threshold=0.0;
      double initial_area;

      std::string file_name;
      if ((dir = opendir(dir_path.c_str())) != NULL) {//work on all files inside the directory(in random order) 
        while ((ent = readdir (dir)) != NULL) {
          file_name=ent->d_name;
          
          if (file_name=="." || file_name=="..") {
            continue;
          }
          std::cout<<"File:"<<file_name<<std::endl;
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
          std::chrono::milliseconds cut_off=std::chrono::milliseconds(500*n_points);
          bool success;  
          
          double anneal_tf;//how much time one annealing iteration takes
          if (preprocess) {
            double anneal_percentage=0.6;
            
            long int base_L=1000;
            double initial_area;
            long int dt=0;
            while (dt==0) {
              base_L=base_L*2;
              Polygon_2 pol;
              Simulated_Annealing<Kernel> anneal;
              auto begin=std::chrono::steady_clock::now();
              anneal.solve(pol,v_points,"min","global","convex_hull",base_L,1,initial_area,cut_off);      
              dt=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-begin).count();

            }
            anneal_tf=(double)dt/(double)cut_off.count();//what percentage one iteration takes out of the cut_off
            glob_L=base_L*std::ceil(anneal_percentage/anneal_tf);//how many iterations should be made so that the total runtime of annealing matches the percentage we want


            base_L=1000;
            dt=0;
            while (dt==0) {
              base_L=base_L*2;
              Polygon_2 pol;
              Simulated_Annealing<Kernel> anneal;
              auto begin=std::chrono::steady_clock::now();
              anneal.solve(pol,v_points,"min","local","convex_hull",base_L,1,initial_area,cut_off);      
              dt=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-begin).count();
            }
            anneal_tf=(double)dt/(double)cut_off.count();//what percentage one iteration takes out of the cut_off
            loc_L=base_L*std::ceil(anneal_percentage/anneal_tf);//how many iterations should be made so that the total runtime of annealing matches the percentage we want



            base_L=1000;
            dt=0;
            while (dt==0) {
              base_L=base_L*2;
              Polygon_2 pol;
              Simulated_Annealing<Kernel> anneal;
              auto begin=std::chrono::steady_clock::now();
              anneal.sub_division(pol,v_points,"min","convex_hull",base_L,initial_area,cut_off);     
              dt=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-begin).count();
            }
            anneal_tf=(double)dt/(double)cut_off.count();//what percentage one iteration takes out of the cut_off
            subd_L=base_L*std::ceil(anneal_percentage/anneal_tf);//how many iterations should be made so that the total runtime of annealing matches the percentage we want
          }     
          
                


          //Hull+Subdivision+SimulatedAnnealing(local)
          int alg_n=0;
          double area;
          for (auto criteria:Min_Max) {
            Polygon_2 polygon_alt;
            Simulated_Annealing<Kernel> anneal;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if ((success=anneal.sub_division(polygon_alt,v_points,criteria,"convex_hull",subd_L,initial_area,cut_off))) { 
              area=abs(polygon_alt.area());
            
              auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
              if (std::chrono::system_clock::now()<(start+cut_off)) {
                cut_off-=milliseconds;
                if ((success=anneal.solve(polygon_alt,v_points,criteria,"local","convex_hull",loc_L,1,initial_area,cut_off))) {
                  area=abs(polygon_alt.area());
                }
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
            Polygon_2 polygon_alt;
            Simulated_Annealing<Kernel> anneal;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if ((success=anneal.sub_division(polygon_alt,v_points,criteria,"convex_hull",subd_L,initial_area,cut_off))) { 
              area=abs(polygon_alt.area());
              auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
              if (std::chrono::system_clock::now()<(start+cut_off)) {
                cut_off-=milliseconds;
                LocalSearch<Kernel> loc(polygon_alt,criteria);
                
                if (!loc.MinimizePolygon(search_L, threshold,cut_off,K)) {
                  area = loc.getPolygonArea();
                  
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
            Polygon_2 polygon_alt;
            Simulated_Annealing<Kernel> anneal;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if (anneal.solve(polygon_alt,v_points,criteria,"global","incremental",glob_L,1,initial_area,cut_off)) {
              area=abs(polygon_alt.area());
            
              auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
              if (std::chrono::system_clock::now()<(start+cut_off)) {
                cut_off-=milliseconds;
                LocalSearch<Kernel> loc(polygon_alt,criteria);
                if (!loc.MinimizePolygon(search_L, threshold,cut_off,K)) {
                  area = loc.getPolygonArea();
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

          
          
          //Hull+SimulatedAnnealing(local)+LocalSearch
          alg_n++;
          for (auto criteria:Min_Max) {
            Polygon_2 polygon_alt;
            Simulated_Annealing<Kernel> anneal;
            cut_off=std::chrono::milliseconds(500*n_points);
            auto start=std::chrono::system_clock::now();
            if ((success=anneal.solve(polygon_alt,v_points,criteria,"local","convex_hull",loc_L,1,initial_area,cut_off))) {
              area=abs(polygon_alt.area());
              auto milliseconds=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
              if (std::chrono::system_clock::now()<(start+cut_off)) {
                cut_off-=milliseconds;
                LocalSearch<Kernel> loc(polygon_alt,criteria);
                if (!loc.MinimizePolygon(search_L,threshold,cut_off,K)) {
                  area = loc.getPolygonArea();
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
          std::cout<<"Done with file:"<<file_name<<std::endl;
        }
        write_results(outfile_path);
        closedir(dir);
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
    if(argv[i][0] == '-' && i < argc) {
      std::string arg = argv[i];
      

      if (arg == "-i") {
        std::string opt = argv[i + 1];
        dir_path=opt;
      }
      else if (arg == "-o") {
        std::string opt = argv[i + 1];
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
  rp.evaluate(dir_path,outfile_path,preprocess);
  return 0;
}
#endif