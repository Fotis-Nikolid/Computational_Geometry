#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <list>
#include <regex>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include "hull.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
int main(int argc,char* argv[]) {
    DIR *dir;
    struct dirent *ent;
    
    std::ifstream file;

    int start=atoi(argv[1]);
    int end=atoi(argv[2]);
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
                        break;
                    }
                }
            }
            std::string file_path=full_path+"/"+file_name;
            file.open(file_path);
            std::string line;
            
            getline(file,line);
            getline(file,line);
            double area;
            if (std::regex_search(line,match_obj,std::regex("\"([0-9]*)\""))) {
                area=atoi(match_obj[1].str().c_str());
            }
            if (std::regex_search(line,match_obj,std::regex("\"([0-9]*)\""))) {
                area=atoi(match_obj[1].str().c_str());
            }
        
            std::list<Point_2> l_points;
            std::vector<Point_2> v_points;
            while (getline(file,line)) {
                double x,y;
                int row;
                std::stringstream sline(line);
                sline >> row;
                sline >> x;
                sline >> y;

                v_points.push_back(Point_2(x,y));
                l_points.push_back(Point_2(x,y));
            }

            //Incremental<K> agl1;
            Hull<K> alg2;

            Polygon_2 p1,p2,p3;
            double random,min,max;
            bool success=true;
            
            random=alg2.solve(p1,l_points,'1');
            if (!p1.is_simple()) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Polygon not simple, 1"<<std::endl;
                perror("Polygon not simple, 1");
                exit(EXIT_FAILURE);
            }
            if (p1.vertices().size()!=points) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Polygon has less points "<<p1.vertices().size()<<"<"<<points<<", 1"<<std::endl;
                exit(EXIT_FAILURE);
            }

            min=alg2.solve(p2,l_points,'2');
            if (!p2.is_simple()) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Polygon not simple, 2"<<std::endl;
                exit(EXIT_FAILURE);
            }
            if (p2.vertices().size()!=points) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Polygon has less points "<<p2.vertices().size()<<"<"<<points<<", 2"<<std::endl;
                exit(EXIT_FAILURE);
            }
            max=alg2.solve(p3,l_points,'3');
            if (!p3.is_simple()) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Polygon not simple, 3"<<std::endl;
                exit(EXIT_FAILURE);
            }
            if (p3.vertices().size()!=points) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Polygon has less points "<<p3.vertices().size()<<"<"<<points<<", 3"<<std::endl;
                exit(EXIT_FAILURE);
            }
            if (min>max) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Min:"<< min <<" greater than Max:"<< max<<std::endl;
                exit(EXIT_FAILURE);
            }
            if (max>area) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" Max:"<< max <<" greater than Area:"<< area<<std::endl;
                exit(EXIT_FAILURE);
            }
            if (random>max || random<min) {
                success=false;
                std::cout<<"Hull error: "<<file_name<<" "<<random <<" random out of range "<<std::endl;
                exit(EXIT_FAILURE);
            }
            if (success) {
                std::cout<<file_name<<" (Hull) ----- OK"<<std::endl;
            }

            

        }
        closedir (dir);
    } 
    else {        // could not open directory 
        perror ("");
        return EXIT_FAILURE;
    }
}

