from PIL import Image, ImageDraw
import os
import json

def print_progress(progress, length, description):
    print("\r%s: %d %% " % (description, int(100*progress/length)), end = "")

def print_progress_unknown(description):
    print("\r%s ... " % (description), end = "")

def print_progress_done(description):
    print("\r%s: Finished " % (description), end = "\n")

def EDGE(path):



    PATH = path

    image = Image.open(PATH)
    SIZE = image.size

    X_ROTATION_AXIS = SIZE[0]/2

    im_rgb = image.convert("RGB")

    #HAS TO BE MODIFIED BEFORE USING REAL IMAGES
    def get_ball(im,width,height):
        image_matrix, ball_points = [], []
        for x in range(width):
            image_matrix.append([]) #Create a new column
            for y in range(height):
                pixel = im.getpixel((x,y))

                #A way of determing if a pixel is ball or non-ball suitable for real photos has to be created here.
                #The four lines below are essentially useless with real photos
                pixel_is_boll = True
                for coulor in pixel:
                    if coulor > 127:
                        pixel_is_boll = False

                image_matrix[x].append(pixel_is_boll) #Append Bool to the newly created column
                if pixel_is_boll:
                    ball_points.append([x,y])

        return image_matrix, ball_points


    def get_mid(ball_points):
        mid_x = 0
        mid_y = 0
        for i in ball_points:
            mid_x += i[0]
            mid_y += i[1]
        mid_x /= len(ball_points)
        mid_y /= len(ball_points)
        return [mid_x, mid_y]

    def is_edge(point,matrix):
        is_edge = False

        #Here we check a square of 9 pixels surrouning a point. If one or more of them is NOT ball then our point is an edge point.
        for x in range(3):
            x_to_check = point[0]-1+x
            if x_to_check > 0 and x_to_check < SIZE[0]: #We make sure we are inside the frame
                for y in range(3):
                    y_to_check = point[1]-1+y
                    if y_to_check > 0 and y_to_check < SIZE[1]: #We maake sure we are inside the frame
                        if not matrix[x_to_check][y_to_check]:
                            is_edge = True
                            return is_edge
                    else:
                        is_edge = True
                        return is_edge
            else:
                is_edge = True
                return is_edge


    def find_edge_layers(ball_points, matrix):
        edge_layers = {}
        for point in ball_points:
            if is_edge(point, matrix):
                if point[1] in edge_layers:
                    edge_layers[point[1]].append(point[0])
                else:
                    edge_layers[point[1]] = [point[0]]

        for l in edge_layers:
            l_sort = sorted(edge_layers[l])
            edge_layers[l] = [l_sort[0]-X_ROTATION_AXIS, l_sort[-1]-X_ROTATION_AXIS]

        return edge_layers

    image_matrix, ball_points = get_ball(im_rgb, SIZE[0], SIZE[1])



    #midpoint = get_mid(ball_points)
    edge_layers = find_edge_layers(ball_points, image_matrix)
    edge_layers_sorted = dict(sorted(edge_layers.items(), key=lambda item: item[0]))

    return edge_layers_sorted

def find_ALL_IMAGES(image_directory):
    images = []
    with os.scandir(image_directory) as dirs:
        dirs = sorted( map( lambda a: a.name, dirs ), key = lambda a : float(a.split(".")[0]) )
        for entry in dirs:
            if entry != '.DS_Store':
                images.append(entry)

    ALL_IMAGES = []

    task_name = 'Scannning images'
    for n, f in enumerate(images):
        print_progress(n, len(images), task_name)
        edge = EDGE(image_directory+images[n])
        ALL_IMAGES.append(edge)

    return ALL_IMAGES

def edge(source, target):

    # Source = filepath to the processed images
    # target = filepath + filename.json where the json should be saved


    image_dir = source
    output_dir = target

    print("\nFIND EDGE")

    write = json.dumps(find_ALL_IMAGES(image_dir))
    with open(output_dir, "w") as out:
        out.write(write)
