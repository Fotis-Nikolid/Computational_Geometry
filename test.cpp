

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
    if ((dir = opendir(path.c_str())) != NULL) {    /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            file_nam=ent->d_name;
            std::cmatch match_obj;
            if (path=="uniform") {
                if (std::regex_search(file.c_str(),match_obj,std::regex("uniform-[0]*([1-9][0-9]*-.*)"))) {
                    points=atoi(match_obj[1]);

                    if (points<start ) {
                        continue
                    }
                    if (points>end) {
                        break;
                    }
                }
            }  
            else {
                if (std::regex_search(file.c_str(),match_obj,std::regex(".*-[0]*([1-9][0-9]*).instance"))) {
                    points=atoi(match_obj[1]);

                    if (points<start ) {
                        continue
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
            if (std::regex_search(line.c_str(),match_obj,std::regex("# parameters \"convex_hull\": {\"area\": \"([0-9]*)\"}"))) {
                area=atoi(match_obj[1]);
            }
        
            std::list<Point_2> l_points;
            std::vector<Point_2> v_points;
            while (getline(file,line)) {
                double x,y;
                int row;
                std::stringstream sline(line);
                sline >> row
                sline >> x;
                sline >> y;

                v_points.push_back(Point_2(x,y));
                l_points.push_back(Point_2(x,y));
            }

            //Incremental<CGAL::Exact_predicates_inexact_constructions_kernel> agl1;
            Hull<CGAL::Exact_predicates_inexact_constructions_kernel> alg2;

            Polygon_2 poly;
            double random,min,max;
            bool success=true;
            
            random=alg2.solve(poly,l_points,"1");
            if (!poly.simple()) {
                success=false;
                std::cout<<"File error: "<<file_name<<std::endl;
                perror("Polygon not simple, 1");
                exit(EXIT_FAILURE);
            }

            min=alg2.solve(poly,l_points,"2");
            if (!poly.simple()) {
                success=false;
                std::cout<<"File error: "<<file_name<<std::endl;
                perror("Polygon not simple, 2");
                exit(EXIT_FAILURE);
            }
            max=alg2.solve(poly,l_points,"3");
            if (!poly.simple()) {
                success=false;
                std::cout<<"File error: "<<file_name<<std::endl;
                perror("Polygon not simple, 3");
                exit(EXIT_FAILURE);
            }
            if (min>max) {
                success=false;
                std::cout<<"File error: "<<file_name<<" "<< min <<" greater than "<< max<<std::endl;
                perror("Min greater than Max");
                exit(EXIT_FAILURE);
            }
            if (max>area) {
                success=false;
                std::cout<<"File error: "<<file_name<<" "<< max <<" greater than "<< area<<std::endl;
                perror("Max greater than Upper Limit");
                exit(EXIT_FAILURE);
            }
            if (random>max || random<min) {
                success=false;
                std::cout<<"File error: "<<file_name<<" "<<random <<" radom out of range "<<std::endl;
                perror("Min greater than Max");
                exit(EXIT_FAILURE);
            }
            if (success) {
                std::cout<<file_name<<"       OK"<<std::endl;
            }

            

        }
        closedir (dir);
    } 
    else {        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
    }
}