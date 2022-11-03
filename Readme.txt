Incremental :

    - First we sort the points with the c++ std::sort function
        - 4 compare function were created , depending if we want to sort for x or y
        - In the case we sort for x , if the points have equal x then we check the y that ensures that we wont
          have points intersect with the red edges
        - The same happens when we sort for y (if they are equal we check the x)
    - Then we create the Triangle with the Initialize function
        - we pop 3 elements from the back of the vector and add them to poylgon
        - if the 3 points have the same x or y then we add more points at the end of the polygon until teh last point 
          have different x and y from the rest (if we did not do that then we would not have a polygon we would have a line)
        - Lastly we set the Convex Hull Polygon the same as the Polygon we have created
    - Then we apply the algorithm bellow until the vector with the points becomes empty
    - First we find where the red edges start and end and we create the new Convex Hull Polygon
        - the red edges are always one after the other ,so if we find some and then we find one edge that is not a red
          edge that means no other red edges exist in the Convex Hull
        - we find where the red edges start and end and we save it in vertices.first_vertex and vertices.second_vertex
        - to find the red edges we use the red_visible function that has the same complexity as the orientation function
          from cgal (the for loop inside the function will loop 3 times max)
        - we add the new vertex in the convex hull polygon between the red edges start and end 
          ( edges : (start , new_vertex) and (new_vertex , end)) And we remove the verteces that are between the start
          and the end (excluding the new one)
        - we return the object vertices (verteces.first_vertex , vertices.second_vertex)
    - Using the object we return obove , we find the candinate edges to be removed from the polygon , we check their
      visibility to new point and we choose one based on the criteria to be removed
        - the way the convex hull polygon and the polygon are constructed , we ensure that if we have 2 verteces
          A and B , and we find A before B in the convex hul polygon ,then if the verteces exist in the polygon 
          we will find A before B there too
        - so we find the verteces.first_vertex and vertices.second_vertex in the polygon and the edges between them
          are the candinate edges to be broken for the new vertex to be added
        - we check if they are visible and based on the criteria (random , min area , max area) we choose 
          one of the edges , and we insert the new vertex between the vertices that create the edge we found