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
#include "polygon.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
int main(int argc,char* argv[]) {//receives starting file number, ending file number and a directory, and performs selected algorithm with all criteria for all files
    DIR *dir;
    struct dirent *ent;
    
    

    int start=atoi(argv[1]);
    int end=atoi(argv[2]);
    std::string algorithm(argv[4]);
    std::string initialization(argv[5]);
    std::string step_choice;
    if (argc==7) step_choice=argv[6];
    std::string path=argv[3];

    std::string full_path="data/"+path;
    int points;

    std::string file_name;
    if ((dir = opendir(full_path.c_str())) != NULL) {    // print all the files and directories within directory 
        while ((ent = readdir (dir)) != NULL) {
            file_name=ent->d_name;
            if (file_name=="." || file_name=="..") {
                continue;
            }
            std::smatch match_obj;
            if (path=="uniform") {
                if (std::regex_search(file_name,match_obj,std::regex("uniform-[0]*([1-9][0-9]*)-[1-9].instance"))) {
                    points=atoi(match_obj[1].str().c_str());

                    if (points<start ) {
                        continue;
                    }
                    if (points>end) {
                        continue;
                    }
                }
            }  
            else {
                if (std::regex_search(file_name,match_obj,std::regex(".*-[0]*([1-9][0-9]*).instance"))) {
                    points=atoi(match_obj[1].str().c_str());

                    if (points<start ) {
                        continue;
                    }
                    if (points>end) {
                        continue;
                    }
                }
            }
            std::string file_path=full_path+"/"+file_name;
            std::ifstream file;
            file.open(file_path);
            if (!file.is_open())
            {
                std::cerr << "There was a problem opening the input file!\n";
                return -1;
            }
            std::string line;
            
            getline(file,line);
            getline(file,line);
            double area;
            if (std::regex_search(line,match_obj,std::regex("\"([0-9]*)\""))) {
                area=atof(match_obj[1].str().c_str());
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


            bool success=true;
            
            Polygon<K> poly(v_points,algorithm,"max",step_choice,initialization,6000,0,10,1);

            if (!poly.Success()) {
                success=false;
                std::cout<<"Polygon Creation Failed..."<<std::endl;
            }
            else {
                Polygon_2 real_poly(poly.get_Polygon());
                if (poly.Size() != v_points.size()) {
                    success=false;
                    std::cout<<"Error: Wrong Number of points "<<poly.Size()<<std::endl;
                }
                if (!poly.Simple()) {
                    success=false;
                    std::cout<<"Error: Not simple"<<std::endl;
                }
                if (poly.Area()<poly.Init_Area()) {
                    success=false;
                    std::cout<<"Error: Bad Optimization "<<poly.Area()<<" vs "<<poly.Init_Area()<<std::endl;
                }
            }
            if (success) {
                std::cout<<"------OK------"<<file_name<<std::endl;
            }
            else {
                std::cout<<"----ERROR------"<<file_name<<std::endl;
            }
            std::cout<<std::endl;
            
            

        }
        closedir (dir);
    } 
    else {        // could not open directory 
        perror ("");
        return EXIT_FAILURE;
    }
}
#endif