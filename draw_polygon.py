import matplotlib.pyplot as plt
  

with open("steps.txt", "r") as file:
    points_x=[]
    points_y=[]
    fl=False
    turn=0
    flag=False
    points=[]
    for str in file:
        if not flag:
            if str[0]=="-":
                flag=True
                continue
            list=str.split(" ")
                
            pointx=float(list[0])

            pointy=float(list[1][:-1])
            points.append([pointx,pointy])
        if flag:
            if str[0]=="-":
                turn+=1
                fl=True
                print(turn)
                plt.plot(points_x,points_y)
                for point in points:
                    plt.plot(point[0],point[1],"o-")
                plt.show()
                points_x.clear()
                points_y.clear()
            else:
                fl=False
                list=str.split(" ")
                
                point1x=list[0]

                point1y=list[1]
                point2x=list[2]
                point2y=list[3][:-1]

                point1=[point1x,point1y]
                point2=[point2x,point2y]
                points_x.append(float(point1[0]))
                points_y.append(float(point1[1]))
                
                points_x.append(float(point2[0]))
                points_y.append(float(point2[1]))
    if (not fl):
        plt.plot(points_x,points_y)
        plt.show()
        print(points_x)
            