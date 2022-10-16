import matplotlib.pyplot as plt
  

with open("steps.txt", "r") as file:
    points_x=[]
    points_y=[]
    fl=False
    for str in file:
        if str[0]=="-":
            fl=True
            plt.plot(points_x,points_y)
            
            plt.show()
            points_x.clear()
            points_y.clear()
            
        else:
            fl=False
            list=str.split(" ")
            point1=list[0]
            point2=list[1]

            point1=point1.split(",")
            point2=point2.split(",")
            points_x.append(float(point1[0]))
            points_y.append(float(point1[1]))
            
            points_x.append(float(point2[0]))
            points_y.append(float(point2[1]))
    if (not fl):
        plt.plot(points_x,points_y)
        plt.show()
        print(points_x)
            