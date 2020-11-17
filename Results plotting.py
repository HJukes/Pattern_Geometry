import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


filenames = ["7 buoy chevron passive tacgem.txt", "Standard bucket 6 buoy passive tacgem.txt"]

for filename in filenames:
    
    data = np.loadtxt(filename, delimiter=',')
    largest = [0,0,0]
    xs = []
    ys = []
    greens = []
    oranges = []
    unique_xs = [data[0][0]]
    greens_for_averaging = [0]
    oranges_for_averaging = [0]
    largest_greens = []
    largest_oranges = []
    counter = 0
    largest_green = 0
    largest_orange = 0
    for i in data:
    # if (i[2] > largest[2]):
    #     largest = i
        
        if (i[0] == unique_xs[counter]) :
        #If the current x is already in the unique xs list, add the z value to the list
            greens_for_averaging[counter] += i[2]
            oranges_for_averaging[counter] += i[3]
            if (i[2] > largest_green):
                
                if (type(i[2]) == list):
                    
                    largest_green = i[2][0]
                else:
                    largest_green = i[2]
            if (i[3] > largest_orange):
                if(type(i[3]) == list):
                    largest_orange = i[3][0]
                else:
                    largest_orange = i[3]
        else:
        #Otherwise if this is the first instance of a new x
            counter +=1
            unique_xs.append(i[0])
            greens_for_averaging.append(i[2])
            oranges_for_averaging.append(i[3])
            if (type(largest_green) == list):
                largest_greens.append(largest_green[0])
            else:
                largest_greens.append(largest_green)
            if (type(largest_orange) == list):
                largest_oranges.append(largest_orange[0])
            else:
                largest_oranges.append(largest_orange)
            
            largest_orange = [i[3]]
            largest_green = [i[2]]
        xs.append(i[0])
        ys.append(i[1])
        greens.append(i[2])
        oranges.append(i[3])
    if (type(largest_green) == list):
        largest_greens.append(largest_green[0])
    else:
        largest_greens.append(largest_green)
    if (type(largest_orange) == list):
        largest_oranges.append(largest_orange[0])
    else:
        largest_oranges.append(largest_orange)
    averaged_greens = []
    averaged_oranges = []
    max_orange = []
    for i in greens_for_averaging:
        averaged_greens.append(i/18)
    for i in oranges_for_averaging:
        averaged_oranges.append(i/18)
    
    

    fig=plt.figure()
    plt.ylabel("Average Kill Zone size")
    plt.xlabel("Angle (°)")
    plt.title("Angle vs Average Kill Zone size for " + filename)
    plt.plot(xs, greens, label="All values")
    plt.plot(unique_xs, averaged_greens, label="Average values", color='green')
    plt.legend()

    fig=plt.figure()
    plt.ylabel("Average Kill Zone size")
    plt.xlabel("Angle (°)")
    plt.title("Angle vs Average orange size for " + filename)
    plt.plot(xs, oranges, label="All values")
    plt.plot(unique_xs, averaged_oranges, label="Average values")
    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.set_ylabel("Orange size")
    ax1.plot(unique_xs, largest_oranges, label="Orange values", color='orange')
    plt.xlabel("Angle (°)")
    plt.title("Angle vs largest size for " + filename)
    #plt.legend()
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("Green size")
    ax2.plot(unique_xs, largest_greens, label="Green values", color='green')
    #plt.legend()

# fig=plt.figure()
# plt.ylabel("Green area",fontsize=16,labelpad=25)
# plt.xlabel("Spacing (x detection range)", fontsize=16)
# plt.title("Spacing vs Kill Zone size for bucket (2 top) 6 buoys")
# plt.plot(ys, zs)

