Incremental :

    - Firstly we sort the points by using c++ std::sort function
        - 4 compare functions have been created for use, depending on whether we want to sort for x or y in descending or ascending order
        - In the case that two points have the same coordinate used for sorting, we also sort based on the other coordinate as to maintain a proper order when inserting           points(otherwise a case occurs that would put a point exactly on top of the existin convex hull, and thus produce 0 visible edges 
    - Afterwards, we create the Triangle with the Initialize function
        - we pop 3 elements from the back of the vector and add them to poylgon(popping from the back is to ensure faster complexity for removing elements from vector)
        - if the 3 points have the same x or y then we add more points at the end of the polygon until the last point 
          has a different x or y from the rest (if we did not do that then we would not have a triangle, we would have a line)
        - Lastly we set the Convex Hull Polygon the same as the Polygon we have created
    - Then we apply the algorithm bellow until the vector with the points becomes empty
    - First we find where the red edges start and end, which we use to build the new Convex Hull Polygon
        - the red edges are always one after the other,so if we find some and then we find one edge that is not red 
          that means no other red edges exist in the Convex Hull and there is no point
        - we find where the red edges start and end and we save it in vertices.first_vertex and vertices.second_vertex
        - to find the red edges we use the red_visible function that utilizes the orientation function
          from cgal, which figures out the orientation of a point compared to a line by using determinants, as is shown in the class lectures(there exists a for loop in           the function, but it will at most be run 3 times, as to find a point in the polygon that is not on the edge that we check for visibility)
        - the red_visible function works in O(1) complexity, and as such, encapsulates the meaning of the incremental algorithm, which allows for incredible fast                   narrowing down of the candidate edge by utilizing the fast lookup for red edges. 
        - we then remove the red edges from the convex hull polygon, and insert the new point between the removed edges
        - the reason we create the new convex hull manually(other than that it is simple to do), is so that we maintain the same orientation between the convex hull               polygon and the real polygon inside of it(this is required to find the edges of the inner polygon)
        - we return the object vertices (verteces.first_vertex , vertices.second_vertex)
    - Using the object we return obove(the two violet points as mentioned in the algorithm presentation), we find the candidate edges to be removed from the polygon 
        - thanks to the way the convex hull polygon and the polygon are constructed , we ensure that if we have 2 vertices
          A and B , and find A before B in the convex hul polygon, then if the verteces exist in the polygon 
          we will find A before B there too.
        - so we find the verteces.first_vertex and vertices.second_vertex in the polygon and the edges between them
          are the candidate edges, for one to be broken for the new point to be added
        - we check if they are visible and based on the criteria (random , min area , max area) we choose 
          one of the edges , and we insert the new vertex between the vertices that create the edge we found
        - the check for an edge visibilty's to a point is O(2k) where k is the number of edges in the real polygon. 
