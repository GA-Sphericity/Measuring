# blender --background --python script.py

import os
import os.path
from os.path import isfile, join
import sys
import bpy, bmesh, os
from math import radians, sin, cos, sqrt, pi
import json
from time import perf_counter, gmtime, strftime, mktime, time

start_path = "" # the starting point for searching files
OB = bpy.data.objects

#### VECTOR MATH #####

def normalize(vector):
    x, y = vector[0], vector[1]
    scale = sqrt(1/(x**2+y**2))
    x *= scale
    y *= scale
    return [x, y]


def rotate_vector(vector, angle):
    sina = sin(angle)
    cosa = cos(angle)
    x = cosa*vector[0] - sina*vector[1]
    y = sina*vector[0] + cosa*vector[1]
    return [x, y]

def rotate_90deg(vector):
    x = -vector[1]
    y = vector[0]
    return [x, y]

def rotate_line(line, angle):
    new_line = [rotate_vector(line[0], angle),rotate_vector(line[1], angle)]
    return new_line


def get_point(line, t):
    px, py, vx, vy = line[0][0], line[0][1], line[1][0], line[1][1]

    x = px + vx * t
    y = py + vy * t
    return [x, y]


def intersection(line0, line1):
    px0, py0, vx0, vy0 = line0[0][0], line0[0][1], line0[1][0], line0[1][1]
    px1, py1, vx1, vy1 = line1[0][0], line1[0][1], line1[1][0], line1[1][1]

    ''' from line0 we get:
    x = px0 + t0 * vx0
    y = py0 + t0 * vy0
    and from line1 we get:
    x = px1 + t1 * vx1
    y = py1 + t1 * vy1

    Solving the equation system for t0 gives us: '''

    t0 = (px1*vy1 - px0*vy1 - vx1*py1 + vx1*py0) / (vx0*vy1 - vy0*vx1)

    '''NB: if our lines are parallel, the divisor = 0 and we will get an error.
    This function should obviously not be used with parallel lines.'''

    return get_point(line0, t0)

def get_line_between_points(point1, point2):
    vector = [point2[0]-point1[0], point2[1]-point1[1]]
    vector = normalize(vector)
    line = [point1, vector]
    return line

def distance_points(point1, point2):
    dist =  sqrt( (point2[0]-point1[0])**2 + (point2[1]-point1[1])**2 )
    return dist

def distance_point_to_line(point, line):
    '''We construct a new line, perpendicular to our
    original line, through our orignal point. We then
    calculate the distance between our originl point and
    the intersection between our new and original lines.'''

    new_line_vector = rotate_90deg(line[1])
    new_line = [point, new_line_vector]

    new_point = intersection(line, new_line)
    dist = distance_points(point, new_point)
    return dist

def is_between_lines(point, line1, line2):

    '''We only use this function when the lines are parallel to one another.
    This means that the distance between the lines can be assumed
    to be constant'''

    distance_between_lines = distance_point_to_line(line1[0], line2)

    dist_to_line1 = distance_point_to_line(point, line1)
    dist_to_line2 = distance_point_to_line(point, line2)

    inside = ( (dist_to_line1 <= distance_between_lines) and (dist_to_line2 <= distance_between_lines) )
    '''If our point is between the lines, the distance
    from the point to each line has to be shorter than
    the distance between the lines'''

    shorter_to_line1 = (dist_to_line1 <= dist_to_line2)
    '''We check which line our point is closest to. True if line1, False if line2'''

    return {'inside':inside, 'shorter_to_line1' : shorter_to_line1}

def get_indices(obj, vertex_group_name):

    vg_idx = obj.vertex_groups[vertex_group_name].index

    vs = [ n for n, v in enumerate(obj.data.vertices) if vg_idx in [ vg.group for vg in v.groups ] ]
    es = []
    for e in obj.data.edges:
        truth_table = []
        for v in e.vertices:
            truth_table.append(vg_idx in [ vg.group for vg in obj.data.vertices[v].groups ])
        if not False in truth_table:
            es.append([v for v in e.vertices])

    return vs, es


#### BOOL2D and BOOL3D : MESH CREATION #####

def bool_2d(obj, vertices_indices, edges_indices, x1, x2, angle, z):

    #we define our lines
    line1 = [rotate_vector([x1,0], angle), rotate_vector([0,1], angle)]
    line2 = [rotate_vector([x2,0], angle), rotate_vector([0,1], angle)]

    new_vertices_index = [] #Vertices to keep is stored here, as indexes
    new_vertices_coordinates = [] #Vertices to create is stored here, as xy coordinates

    verts = obj.data.vertices

    for e in edges_indices:

        line_relations_list = []
        vert_co = [] #Here we store the xy coordinates of each vertex of an edge

        for v in e:

            p = [ verts[v].co[0], verts[v].co[1] ] #we get xy coords
            vert_co.append(p)

            if not (v in new_vertices_index): #If we havent already looked at this vertice and found that it is inside
                line_relations = is_between_lines(p, line1, line2)
                if line_relations['inside']:
                    new_vertices_index.append(v)

                line_relations_list.append(line_relations)

            else:
                line_relations_list.append({'inside':True})



        if (not line_relations_list[0]['inside']) and (not line_relations_list[1]['inside']):

            if (line_relations_list[0]['shorter_to_line1'] != line_relations_list[1]['shorter_to_line1']):

                new_line = get_line_between_points(vert_co[0], vert_co[1])
                inter0 = intersection(line1, new_line)
                inter1 = intersection(line2, new_line)

                new_vertices_coordinates.append(inter0)
                new_vertices_coordinates.append(inter1)

        elif (line_relations_list[0]['inside'] != line_relations_list[1]['inside']):

            outside_vert = int(line_relations_list[0]['inside'])

            if line_relations_list[outside_vert]['shorter_to_line1']:
                line_to_use = line1
            else:
                line_to_use = line2

            #Construct a linear function between the two known vertices
            new_line = get_line_between_points(vert_co[0], vert_co[1])

            #Find the intersection of our original function and our new function
            inter = intersection(line_to_use, new_line)

            new_vertices_coordinates.append(inter)

    for v in new_vertices_index:
        p = [ verts[v].co[0], verts[v].co[1] ]
        new_vertices_coordinates.append(p)


    for v in new_vertices_coordinates:
        v.append(z)

    return new_vertices_coordinates

def bool_3d(obj, ob2, layer_lines, angle):

    bpy.ops.object.mode_set(mode='OBJECT')

    layers = {}

    #### CALCULATE OUR NEW LAYERS ####

    for g in obj.vertex_groups:
        layer_key = g.name
        if layer_key in layer_lines:
            verts, edges = get_indices(obj, g.name)
            if len(verts)  > 0:
                z = obj.data.vertices[verts[0]].co[2]
                x_lines = layer_lines[layer_key]
                layer_co = bool_2d(obj, verts, edges, x_lines[0], x_lines[1], angle, z)
                layers[layer_key] = layer_co



    #### WE WILL NOW MANIPULATE OUR MESH ####

    create_object(layers, ob2)

def create_object(layers, ob2):

    layer_indices = {}

    me = ob2.data

    bm = bmesh.new()
    bm.from_mesh(me)
    bm.clear()

    current_vert = 0
    for l in layers:
        #print("creeate_object: layers: ", layers[l])
        layer_indices[l] = []
        for v in layers[l]:
            #print("v: ", v)
            bm.verts.new(v)
            layer_indices[l].append(current_vert)
            current_vert +=1

    bm.to_mesh(me)
    bm.free

    #### REMOVING OLD VERTEX GROUPS AND ADDING NEW ONES ####

    ob2.vertex_groups.clear()

    for l in layer_indices:
        group = ob2.vertex_groups.new( name = str(l) )
        group.add( [v for v in layer_indices[l]], 1.0, 'ADD' )

    #### UTILISING BPY.OPS TO ADD FACES (AND THUS EDGES) TO OUR LAYERS ####

    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = ob2
    ob2.select_set(True)

    bpy.ops.object.mode_set(mode='EDIT')

    for l in layer_indices:
        bpy.ops.mesh.select_all(action='DESELECT')
        for v in layer_indices[l]:
            bpy.ops.object.vertex_group_set_active(group=str(l))
            bpy.ops.object.vertex_group_select()

            bpy.ops.mesh.edge_face_add()

    bpy.ops.object.mode_set(mode='OBJECT')


#### IMAGE INETRPRETATION ###

#höhö

def print_progress(progress, length, description):
    print("\r%s: %d %% " % (description, int(100*progress/length)), end = "")

def print_progress_unknown(description):
    print("\r%s ... " % (description), end = "")

def print_progress_done(description):
    print("\r%s: Finished " % (description), end = "\n")



#### MODEL BALL ####

def create_initial_object(all_images, ob2):
    init_im = [0, int(len(all_images)/2)] # HE
    init_angles = [ 0, radians(180*init_im[1]/len(all_images)) ]

    layers = {}


    for l in all_images[ init_im[0] ]: #l is the name of the layer
        if l in all_images[ init_im[1] ]:

            try:
                layer1 = all_images[ init_im[0] ][l]
            except:

                print('l: ', l)
            layer2 = all_images[ init_im[1] ][l]

            line1_1 = [rotate_vector([layer1[0],0], init_angles[0]), rotate_vector([0,1], init_angles[0])]
            line1_2 = [rotate_vector([layer1[1],0], init_angles[0]), rotate_vector([0,1], init_angles[0])]

            line2_1 = [rotate_vector([layer2[0],0], init_angles[1]), rotate_vector([0,1], init_angles[1])]
            line2_2 = [rotate_vector([layer2[1],0], init_angles[1]), rotate_vector([0,1], init_angles[1])]



            p1 = intersection(line1_1, line2_1)
            p2 = intersection(line1_1, line2_2)
            p3 = intersection(line1_2, line2_1)
            p4 = intersection(line1_2, line2_2)

            point_list = [p1, p2, p3, p4]
            point_list_fixed = []

            for p in point_list:
                p.append(float(l))
                point_list_fixed.append(p)

            layers[l] = point_list

    create_object(layers, ob2)
    return init_im



### CALCULATE SPHERICITY ###

def get_sphericity(obj):

    me = obj.data
    bm = bmesh.new()
    bm.from_mesh(me)

    A = sum(f.calc_area() for f in bm.faces)
    V = bm.calc_volume()

    def sphericity(A, V):
        return ((pi**(1/3)) * ( (6*V)**(2/3) ))/A

    return sphericity(A,V)



### RUNNER ###

def create_3d_model(ALL_IMAGES, goal_object):

    bpy.context.view_layer.objects.active = goal_object
    bpy.ops.object.mode_set(mode='OBJECT')

    angle_change = radians(180/len(ALL_IMAGES))
    already_done = create_initial_object(ALL_IMAGES, goal_object)

    task_name = 'Modeling object'
    for i in range(len(ALL_IMAGES) ):
        print_progress(i, len(ALL_IMAGES), task_name)
        if not (i in already_done):
            bool_3d(goal_object, goal_object, ALL_IMAGES[i], i*angle_change)




def analyze(name, ALL_IMAGES):

    print("Analyzing %s in directory:\n%s" % (name, "dir") )

    start_time = perf_counter()

    '''First off, we set up our objects and collections'''

    parent_collection_name = 'Collection'

    no_mod_name = name + '_no_mod'

    mesh_no_mod = bpy.data.meshes.new(no_mod_name)
    object_no_mod = OB.new(no_mod_name, mesh_no_mod)

    mesh = bpy.data.meshes.new(name)
    object = OB.new(name, mesh)

    collection_parent = bpy.data.collections.get(parent_collection_name)

    collection_child = bpy.data.collections.new(name)
    collection_parent.children.link(collection_child)

    collection_child.objects.link(object_no_mod)
    collection_child.objects.link(object)

    '''Then, we add modifiers with the right settings to our main object'''
    solidify = object.modifiers.new(name="Solidify", type='SOLIDIFY')
    solidify.thickness = 1
    solidify.offset = 0

    remesh = object.modifiers.new(name="Remesh", type='REMESH')
    remesh.mode = 'VOXEL'
    remesh.voxel_size = 3

    smooth = object.modifiers.new(name="Smooth", type='SMOOTH')
    smooth.factor = 1
    smooth.iterations = 4


    ''' We create our raw model on our no_mod object '''
    create_3d_model(ALL_IMAGES, object_no_mod)

    ''' We copy mesh data of our object_no_mod to our object with modifiers '''
    task_name = 'Applying modifiers'
    print_progress_unknown(task_name)
    object.data = object_no_mod.data.copy()
    object_no_mod.hide_viewport = object_no_mod.hide_render = True

    bpy.context.view_layer.objects.active = object
    bpy.ops.object.modifier_apply(modifier = solidify.name)
    bpy.ops.object.modifier_apply(modifier = remesh.name)
    bpy.ops.object.modifier_apply(modifier = smooth.name)

    print_progress_done(task_name)

    ''' We calculate the sphericity '''
    task_name = 'Calculating sphericity'
    print_progress_unknown(task_name)

    sphericity = get_sphericity(object)
    print_progress_done(task_name)

    '''We print results '''
    end_time = perf_counter()
    print('Finished in %s' % ( strftime("%H h %M min %S sec", gmtime(end_time - start_time)) ) )
    object["sphericity"] = sphericity
    print('Sphericity: ',sphericity,'\n')
    return sphericity

#### Runner ###
"""
ALL_IMAGES_json = '/home/ingvar/Documents/GA/data/Ball_Json/n2_m1_nr1_analyzed.json'

name = ALL_IMAGES_json.split('/')[-1].split('.')[0]

with open(ALL_IMAGES_json, "r") as aid:
    ALL_IMAGES = json.loads(aid.read())

#print(ALL_IMAGES)

analyze(name, ALL_IMAGES)
"""

# bpy.ops.wm.save_as_mainfile(filepath=filepath)

def remove_obj_and_mesh():

    for collection in bpy.data.collections:
        if 'Data' in collection.name:
            current_collection = bpy.data.collections.get(collection.name)
            bpy.data.collections.remove(current_collection)


starttime = time()

roundness = list()

# Look though every file
for path, directories, files in os.walk(start_path):
    for directory in directories:
        current_path = path + '/' +  directory
        if "Images" == current_path.split('/')[-1]: # Filter files so only image directories show up
            print("=" * 80)
            print(current_path)

            # OBS LOOK IF THERE ALREADY IS A .BLEND FILE, SKIP THAT DIRECTORY IN THAT CASE.
            # OBS LOOK IF THERE IS A EDGES.JSON AVAILABLE

            # Set the save file for the .blend file visulising the ball
            save_data_file = "/".join(current_path.split('/')[:len(current_path.split('/')) - 1]) + "/visual.blend"
            # Set the save file for the .json file containing information about
            save_json_file = "/".join(current_path.split('/')[:len(current_path.split('/')) - 1]) + "/data.json"
            # The path to the edges json file
            edges_file = current_path + "/edges.json"
            print(".blend save dir   :", save_data_file)
            print("json data save dir:", save_json_file)


            if os.path.exists(save_data_file) == False and os.path.exists(edges_file) == True:

                # Gather edges data
                with open(edges_file, "r") as current_edges_file:
                    edges_data = json.loads(current_edges_file.read())

                # The file analysing, the processing heavy operation.
                sphericity = analyze("Data", edges_data)

                roundness.append(sphericity)
                print(roundness)

                # The JSON data dump
                data_file_dump = json.dumps({"sphericity" : sphericity})
                print(data_file_dump)

                with open(save_json_file, "w") as write_json:
                    write_json.write(data_file_dump)

                # Saves the .blend file
                bpy.ops.wm.save_as_mainfile(filepath=save_data_file)

            else:
                # Write out the file 'error'
                if os.path.exists(save_data_file) == True:
                    print(".blend file already exists")
                else:
                    print("Missing edges.json")

            # Remove data so a duplication is not made in the other file.

            remove_obj_and_mesh()

endtime = time()

print("Total time:", endtime - startime)

'''

DIR = '<Directory>'
m = bpy.data.texts['CreateModel'].as_module()
m.analyze('Name', DIR)

'''
