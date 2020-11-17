import numpy as np
import matplotlib.pyplot as plt
#import csv

xy = -20 + 50 * np.random.uniform(size=(10000, 2))

def make_chevron(angle, spacing, radius):
    """
    Makes a list of 5 circles in a Chevron
    """
    cos_angle = np.cos(angle * np.pi/180)
    sin_angle = np.sin(angle * np.pi/180)
    circles = [ # In the format x, y, radius
        [0,                   0,                   radius], # origin circle
        [spacing,             0,                   radius], # first circle to the right
        [2*spacing,           0,                   radius], # second circle to the right
        [spacing*cos_angle,   spacing*sin_angle,   radius], # first diagonal circle up
        [2*spacing*cos_angle, 2*spacing*sin_angle, radius]  # second diagonal circle
    ]
    return circles

def make_circles(angle, spacing, radius):
    """
    Makes a list of 5 circles in a bent 'L' with
    given angle and spacing between circle centres.
    """
    cos_angle = np.cos(angle * np.pi/180)
    sin_angle = np.sin(angle * np.pi/180)
    circles = [
        [0,                   0,                   radius], # origin circle
        [spacing,             0,                   radius], # first circle to the right
        [2*spacing,           0,                   radius], # second circle to the right
        [spacing*cos_angle,   spacing*sin_angle,   radius], # first diagonal circle up
        [2*spacing*cos_angle, 2*spacing*sin_angle, radius]  # second diagonal circle
    ]
    return circles


def make_circles_H(angle, spacing, radius, number):
    cos_angle = np.cos(angle * np.pi/180)
    sin_angle = np.sin(angle * np.pi/180)
    counter = 0
    circles = []
    if (number % 2 == 1): #If the number is odd
        circles.append([0,0,radius])        
        while (counter < ((number -1)/2)):
            circles.append([(counter+1)*spacing, 0, radius])
            circles.append([(counter+1) * cos_angle * spacing, (counter+1) * spacing * sin_angle, radius])
            counter +=1
    
    return circles

def make_bucket_2top(angle, spacing, radius, number):
    cos_angle = np.cos(angle * np.pi/180)
    sin_angle = np.sin(angle * np.pi/180)
    counter = 1
    circles = []
    circles.append([0,0,radius])
    while (counter < (number - 2)):
        circles.append([counter * spacing, 0, radius])
        counter +=1
    circles.append([-sin_angle * spacing, cos_angle * spacing, radius])
    circles.append([(number -3) * spacing + sin_angle * spacing, cos_angle * spacing, radius])
    return circles

def make_bucket_2bot(angle, spacing, radius, number):
    cos_angle = np.cos(angle * np.pi/180)
    sin_angle = np.sin(angle * np.pi/180)
    counter = 1
    circles = []
    circles.append([0,0,radius])
    circles.append([spacing, 0, radius])
    while (counter < ((number - 2)/2)+1):
        circles.append([counter * -1* sin_angle * spacing, counter * cos_angle * spacing, radius])
        circles.append([spacing + (counter * sin_angle * spacing), counter * cos_angle * spacing, radius])
        counter+=1
    return circles

def count_overlaps_fast(xy, circles):
    """
    For each xy_i in xy matrix, counts overlapping circles.
    ARGUMENTS
        xy: (num_points, 2) matrix of points
        circles: list of [x, y, r] triplets
        
    RETURNS:
        counts: (num_points) vector of integers
                showing number of overlapping
                circles at point xy_i
    """
    circles = np.array(circles) # shape (num_circles, 3)
    
    # reshape centres to (2, 1, num_circles) array
    circle_centres = circles[:, :2].reshape(1, -1, 2)
    circle_centres = np.transpose(circle_centres, (2, 0, 1))
    
    # reshape points to (2, num_points, 1) array
    xy = xy.reshape(1, -1, 2)
    xy = np.transpose(xy, (2, 1, 0))
    
    # auto broadcasts to (2, num_points, num_circles)
    SQ_dist = (xy - circle_centres)**2
    
    # summation over first axis, (num_points, num_circles)
    SQ_dist_matrix = np.sum(SQ_dist, 0)
    
    # squared radius of each circle, reshape to (1, num_circles)
    R2 = circles[:, 2]**2
    R2 = R2.reshape(1, -1)
    
    # in each column, count number of SQ_dist within radius
    counts = np.sum( SQ_dist_matrix <  R2, 1)
    
    return counts

def plot_points_and_circles(circles):
    # plot points, circles
    fig, ax = plt.subplots()
    counts = count_overlaps_fast(xy, circles)
    triple_mask = counts > 2
    ax.scatter(xy[~triple_mask, 0], xy[~triple_mask, 1], alpha=0.1)
    ax.scatter(xy[triple_mask, 0], xy[triple_mask, 1], alpha=0.5)

    # add each circle to the plot
    for c in circles:
        circle_plot_object = plt.Circle(c[:2], c[2], color='r', alpha=0.1)
        ax.add_artist(circle_plot_object)
        
def check_circle_angles(xy, circles, max_angle=180):
    """
    Computes the angles from xy to each circle centre, then
    checks inter-circle angle is below max_angle
    """
    
    # not enough neighbors, type = 0
    if circles.shape[0] < 3:
        return 0
    
    circles_centres = np.array(circles[:, :2])
    
    dx = xy[0] - circles_centres[:, 0]
    dy = xy[1] - circles_centres[:, 1]
    
    # get ther angle to all neighbor centres
    angles = np.sort(np.arctan2(dx, dy))
    angles = np.hstack([angles, angles[0] + 2*np.pi])
    
    # get the angle differences
    delta_angles = np.abs(angles[1:] - angles[:-1])
    
    # largest angle separation is over limit, type = 1
    if delta_angles.max() > max_angle * np.pi / 180:
        return 1
    
    else: # largest angle is within limit, type = 2
        return 2

def check_triple_and_angles(xy, circles, max_angle=180):
    """
    Finds all xy points with 3 overlapping circles. Then filters
    those points to get ones whose circles are angularly well spaced,
    i.e. the inter-circle angles are at most max_angle.
    """
    circles = np.array(circles).reshape(-1, 3)
    
    # for each xy point, find the nearest neighbore circles
    nbrs, counts = get_neighbor_circles(xy, circles)
    
    
    # get the point type, 
    #     0: useless
    #     1: triple coverage but not surrounded
    #     2: triple coverage and surrounded
    xy_type = []
    for xy_i, nbr_i in zip(xy, nbrs):
        type_i = check_circle_angles(xy_i, circles[nbr_i, :])
        xy_type.append(type_i)
        
    xy_type = np.array(xy_type)
    
    return xy_type

def plot_points_types(circles, plot=True):
    #plot points, circles
    
    circles = np.array(circles)
    min_x = (circles[:, 0] - circles[:, 2]).min()
    max_x = (circles[:, 0] + circles[:, 2]).max()
    min_y = (circles[:, 1] - circles[:, 2]).min()
    max_y = (circles[:, 1] + circles[:, 2]).max()
    x_points = min_x + (max_x - min_x) * np.random.uniform(size=(100000, 1))
    y_points = min_y + (max_y - min_y) * np.random.uniform(size=(100000, 1))
    xy = np.hstack([x_points, y_points])
    types = check_triple_and_angles(xy, circles)
    green_area = (np.mean(types==2))*(max_x - min_x)*(max_y-min_y)
    orange_area = (np.mean(types==1))*(max_x - min_x)*(max_y-min_y) + green_area
    
    if (plot != False):
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_title("orange area:"+str(orange_area)+", green area:"+ str(green_area))
        ax.set_aspect("equal")
        #ax.legend((str(orange_area), str(green_area)),('orange area', 'green area'))
        for t in np.unique(types):
            mask = t == types
            alpha = [0.1, 0.1, 1][int(t)]
            ax.scatter(xy[mask, 0], xy[mask, 1], alpha=alpha)
        if np.sum(mask)==0: ax.scatter(0, 0, alpha=0)
        #add each circle to the plot
        for c in circles:
            circle_plot_object = plt.Circle(c[:2], c[2], color='r', alpha=0.1)
            ax.add_artist(circle_plot_object)
    
    return green_area, orange_area

def get_neighbor_circles(xy, circles):
    """
    For each xy_i in xy matrix, counts overlapping circles.
    ARGUMENTS
        xy: (num_points, 2) matrix of points
        circles: list of [x, y, r] triplets
        
    RETURNS:
        counts: (num_points) vector of integers
                showing number of overlapping
                circles at point xy_i
    """
    circles = np.array(circles) # shape (num_circles, 3)
    
    # reshape centres to (2, 1, num_circles) array
    circle_centres = circles[:, :2].reshape(1, -1, 2)
    circle_centres = np.transpose(circle_centres, (2, 0, 1))
    
    # reshape points to (2, num_points, 1) array
    xy = xy.reshape(1, -1, 2)
    xy = np.transpose(xy, (2, 1, 0))
    
    # auto broadcasts to (2, num_points, num_circles)
    SQ_dist = (xy - circle_centres)**2
    
    # summation over first axis, (num_points, num_circles)
    SQ_dist_matrix = np.sum(SQ_dist, 0)
    
    # squared radius of each circle, reshape to (1, num_circles)
    R2 = circles[:, 2]**2
    
    # for each row, get the neighbour circles (if any) with distance
    ngbr_circles = [np.where(sq<R2)[0] for sq in SQ_dist_matrix]
    
    counts = np.array([len(nbrs) for nbrs in ngbr_circles])
    
    return ngbr_circles, counts


def Compute_Kill_Zone(circles):
    """
    Circles is a list of lists of radius x and y points. 
    The result is the orange area is the trip coverage and green is trip coverage + 60 degree crosscut.
    Make_Circles = chevron.
    """
    # plot points, circles
    circles = np.array(circles)
    min_x = (circles[:, 0] - circles[:, 2]).min()
    max_x = (circles[:, 0] + circles[:, 2]).max()
    min_y = (circles[:, 1] - circles[:, 2]).min()
    max_y = (circles[:, 1] + circles[:, 2]).max()
    x_points = min_x + (max_x - min_x) * np.random.uniform(size=(100000, 1))
    y_points = min_y + (max_y - min_y) * np.random.uniform(size=(100000, 1))
    xy = np.hstack([x_points, y_points])
    types = check_triple_and_angles(xy, circles)
    green_area = (np.mean(types==2))*(max_x - min_x)*(max_y-min_y)
    orange_area = (np.mean(types==1))*(max_x - min_x)*(max_y-min_y) + green_area
    return green_area, orange_area


def Kill_Grid(spacing, radius, nx, ny):
    """
    Creates a grid of buoys of x and y in number
    """
    circles = []
    for i in range(nx):
        for j in range(ny):
            new_circle = [i*spacing, j*spacing, radius]
            circles.append(new_circle)
            
    return circles

#results = []
min_test_angle_bucket = 0
max_test_angle_bucket = 31
min_test_angle_chevron = 60
max_test_angle_chevron = 121
test_angle_step = 30
min_spacing = 0.5
max_spacing = 1.1
spacing_step = 0.25



test_angles_chevron = np.arange(min_test_angle_chevron, max_test_angle_chevron, step=test_angle_step)
test_angles_buckets = np.arange(min_test_angle_bucket, max_test_angle_bucket, step=test_angle_step)
test_spacings = np.arange(min_spacing, max_spacing, step=spacing_step)

print(test_angles_chevron)
print(test_spacings)

def RunTest(testingangles, testingspacings, nobuoys, run):
    angles = []
    spacings = []
    greens = []
    oranges = []
    for i in testingangles:
        print("run %i" %run)
        print(i)
        for j in testingspacings:
            print(j)
            if (run == 1) or (run == 2):
                green, orange = plot_points_types(make_bucket_2top(i, j, 1, nobuoys))
            elif (run == 3):
                green, orange = plot_points_types(make_bucket_2bot(i, j, 1, nobuoys), False)
            elif (run == 4) or (run == 5):
                green, orange = plot_points_types(make_circles_H(i, j, 1, nobuoys))
            angles.append(i)
            spacings.append(j)
            greens.append(green)
            oranges.append(orange)

    counter = 0
    if (run == 1):
        tit = "Standard bucket 4 buoy"
    elif (run ==2):
        tit = "Standard bucket 6 buoy passive tacgem"
    elif (run == 3):
        tit = "Elongated bucket 6 buoy"
    elif (run == 4):
        tit = "7 buoy chevron passive tacgem"
    elif (run == 5):
        tit = "7 buoy chevron"
    with open(tit + '.txt', 'w') as writer:
        while (counter < len(angles)):
        #results += str(angles[counter]) & "," & str(spacings[counter]) & "," & str(greens[counter]) & "\n"
            writer.write(str(angles[counter]) + "," + str(spacings[counter]) + "," + str(greens[counter]) + "," + str(oranges[counter]) + "\n")
            counter +=1


#plot_points_types(make_bucket_2bot(0, 1, 1, 6))

#Test 1, bucket 2 top, 4 buoys
#RunTest(test_angles_buckets, test_spacings, 4, 1)
#Test 2, bucket 2 top, 6 buoys
#RunTest(test_angles_buckets, test_spacings, 6, 2)
#Test 3, bucket 2 bottom, 6 buoys
#RunTest(test_angles_buckets, test_spacings, 6, 3)
#Test 4, 5 buoy chevron
RunTest(test_angles_chevron, test_spacings, 7, 4)
#Test 5, 7 buoy chevron
#RunTest(test_angles_chevron, test_spacings, 7, 5)



#plot_points_types(make_circles_H(26, 1, 1, 7), True)

#plot_points_and_circles(make_bucket_2bot(-45, 10, 5, 4))

#print(make_bucket_2top(90, 1, 1, 8))

























#plot_points_types(make_circles_H(37,1,1,7))
    
#print(spacings)

#plot_points_types(make_circles_H(80,1,1,5))

#make_circles(90, 10, 10)
#make_circles(60, 10, 10)

#C1 = make_circles(75, 10, 10)
#C1 + [[5, 5, 10]]
#count_overlap(0, 0, C1)
#count_overlap(5, 5, C1)


# make a load of points in the (-20, 30) x (-20, 30) square
#xy = -20 + 50 * np.random.uniform(size=(10000, 2))
#counts = count_overlaps_fast(xy, C1)

#triple_coverage_area = 50 * 50 * np.mean(counts > 2)
#print("Triple coverage area:", triple_coverage_area)


#plot_points_and_circles(C1)
#plot_points_and_circles(make_circles(90, 7.5, 10))
#plot_points_and_circles(make_circles(60, 7.5, 10))
#plot_points_and_circles(make_circles(60, 7.5, 10))

#spacing= np.linspace(7, 8, 3)
#angles= np.linspace(59, 60, 3)
#areas=[]
# for space in spacing:
#     for angle in angles:
#         area = plot_points_types(make_circles(angle, space, 10))
#         area = [space, angle] + [area[0], area[1]]
#         areas.append(area)
        
# print(areas)


# for res in areas:
#     print (res)
# Kill_Grid (1, 1, 10, 3)

# plot_points_types(Kill_Grid(0.99,1,10,2))
# plot_points_types(Kill_Grid(0.9,1,3,3))


# spacing= np.linspace(4, 12, 9)
# angles= np.linspace(30, 100, 8)
# areas=[]
# for space in spacing:
#     for angle in angles:
#         area = Compute_Kill_Zone(make_circles(angle, space, 10))
#         area = [space, angle] + [area[0], area[1]]
#         areas.append(area)
        
# print(areas)
# print("FINISHED")

