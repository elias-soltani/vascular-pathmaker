"""
Utility functions for zinc files.
"""

import numpy as np
import os


from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldCoordinates,\
    findOrCreateFieldStoredString
from opencmiss.zinc.context import Context
from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field, FieldGroup
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK
from opencmiss.zinc.element import MeshGroup
from opencmiss.utils.zinc.general import ChangeManager
from scaffoldmaker.utils.interpolation import getCubicHermiteArcLength, interpolateCubicHermite



def mirror_zinc_file(zinc_file, plane):
    """
    Mirror the zinc file with given plane
    :param zinc_file:
    :param plane:  [n1,n2,n3, d]  n1x+n2y+n3z=d, where n=[n1,n2,n3] is plane normal vector and d is plane parameter.
    only works for [1, 0, 0, 0] for now
    :return:
    """
    dirname = os.path.dirname(zinc_file)
    basename = os.path.basename(zinc_file)
    output_zinc_mirrored = os.path.join(dirname+r'\output', basename.split('_right')[0]+'_left.exf')

    element_started = False
    counter = -100
    # define_started = False
    with open(zinc_file, 'r') as f, open(output_zinc_mirrored, 'w') as g:
        for line in f:
            if element_started:
                g.write(line)
                continue
            # if 'node2' in line:
            #     define_started = True
            if 'Node:' in line:
                counter = 0

            elif 'mesh' in line:
                element_started = True
            else:
                counter += 1
                if counter == 1:  # Todo modify this for general plane.
                    line = line.strip().split()
                    line = [str(-float(c)) for c in line]
                    line = ' '.join(line)
                    counter = -100
            g.write(line)


def group_all_elements_nodes_of_zinc(zinc_files_directory):
    """
    Modify cmgui com file for each zinc file in the folder and output the zinc file with the whole group in the output
    directory one by one.
    :param zinc_files_directory:
    :return:
    """
    import subprocess
    import glob
    cmgui_exe = r'C:\Users\egha355\Desktop\sparc3\Mahyar cmgui\cmgui.exe'
    input_cmgui = os.path.join(zinc_files_directory, 'convert.cmgui')
    assert os.path.isfile(input_cmgui), 'convert.cmgui file was not found in the directory'
    output_cmgui = os.path.join(zinc_files_directory, 'convert2222222222222222222.cmgui')
    if not os.path.isdir(zinc_files_directory + r'\output'):
        os.mkdir(zinc_files_directory + r'\output')
    for file in glob.glob(zinc_files_directory + r'\*'):
        if os.path.isfile(file) and not os.path.basename(file) in ['convert.cmgui', 'convert2222222222222222222.cmgui']:
            print('Import com file convert2222222222222222222.cmgui ')
            with open(input_cmgui, 'r') as f, open(output_cmgui, 'w') as g:
                for line in f:
                    line = line.replace('file_name', os.path.basename(file).split('.')[0])
                    g.write(line)
            subprocess.call([cmgui_exe])
    os.remove(output_cmgui)


def get_vessel_group_names(field_module):
    """

    :param region:
    :param field_module:
    :return:
    """
    fielditer = field_module.createFielditerator()
    field = fielditer.next()

    groups = []

    while field.isValid():
        field_name = field.getName()
        if '.mesh1d' in field_name:
            groups.append(field_name.split('.mesh1d')[0])
        field = fielditer.next()

    return groups


def get_vessel_lengths(field_module, groups):
    """
    Read the zinc file and for each element, calculate the length. We assume that we have cubic Hermits basis functions.
    :param field_module:
    :return:
    """
    lengths = []
    radius1 = []
    radius2 = []
    points = {}

    cache = field_module.createFieldcache()
    coordinates = field_module.findFieldByName('coordinates').castFiniteElement()
    radius_field = field_module.findFieldByName('radius').castFiniteElement()
    # Get element iterator
    mesh = field_module.findMeshByDimension(1)

    dxi = 1.0  # mm
    for group in groups:
        field_group = field_module.findFieldByName(group).castGroup()
        element_group = field_group.getFieldElementGroup(mesh)
        if element_group.isValid():
            mesh_group = element_group.getMeshGroup()
        else:
            print('extractPathParametersFromRegion: missing group "' + group + '"')
        elemiter = mesh_group.createElementiterator()
        element = elemiter.next()
        while element.isValid():
            xv = []
            # for each element use Gaussian quadrature to calculate the arc length
            eft = element.getElementfieldtemplate(coordinates, 3)
            node1 = element.getNode(eft, 1)
            node2 = element.getNode(eft, 2)
            cache.setNode(node1)
            result, r1 = radius_field.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, v1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            cache.setNode(node2)
            result, r2 = radius_field.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, v2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            arcLength = getCubicHermiteArcLength(v1, d1, v2, d2)
            radius1.append(r1)
            radius2.append(r2)
            lengths.append(arcLength)

            # point cloud
            for i in range(int(arcLength // dxi)):
                xv.append(interpolateCubicHermite(v1, d1, v2, d2, i * dxi / arcLength))
            # import matplotlib.pyplot as plt
            # ax = plt.figure().add_subplot(projection='3d')
            # xv=np.array(xv)
            # ax.plot(xv[:, 0], xv[:, 1], xv[:, 2], 'o')
            # plt.show()
            points[group] = xv

            element = elemiter.next()
            if element.isValid():
                print('Warning: More than one element per group')

    # relationships = find_vessels_relationships(groups)
    relationships = []
    distances = []
    import time
    start = time.time()
    for group in groups:
        min_dist1 = 1000
        min_dist2 = 1000
        # find closest point to each end of the vessel and add it to reationships
        for name in groups:
            if name != group:
                dist1 = find_closest_point_to_vessel(points[group][0], points[name])[0]
                dist2 = find_closest_point_to_vessel(points[group][-1], points[name])[0]
                if dist1 <= min_dist1:
                    min_dist1 = dist1
                    name1 = name
                if dist2 <= min_dist2:
                    min_dist2 = dist2
                    name2 = name
        relationships.append([name1, name2])
        distances.append([str(min_dist1), str(min_dist2)])
    print('Time to find relationships: ', time.time() - start)
    return lengths, radius1, radius2, relationships, distances


def find_closest_point_to_vessel(point, points):
    """
    find the closest point to the vessel
    :param point:
    :param points:
    :return:
    """
    from scipy.spatial import cKDTree
    import numpy as np

    point = np.array(point)
    points = np.array(points)

    tree = cKDTree(points)
    dist, idx = tree.query(point)
    return dist, idx


def _write_vessel_geometry(csv_file, groups, lengths, radius1, radius2, relationships, distances):
    with open(csv_file, 'w') as g:
        for i in range(len(groups)):
            g.write(groups[i].split('.mesh1d')[0] + ',' + str(lengths[i]) + ',' + str(radius1[i]) + ',' + str(radius2[i]) + ',' +
                    ','.join(relationships[i]) + ','+ ','.join(distances[i]) + '\n')


def get_groups_with_one_element(groups):
    """
    Get only numbered elements and exclude the tree.
    Some groups have multiple elements which in the arteries file each element of them are names with numbering.
    :param groups: list of strings
    :return: list of strings
    """
    import re
    groups_excluded = [c for c in groups]
    for name in groups:
        if re.search(r' 1["\s]', name):
            groups_excluded.remove(name.replace(' 1', ''))

    return groups_excluded


def output_vessel_anatomical_properties(scaffold_file):
    """

    :param scaffold_file:
    :return:
    """
    # region = load_zinc_file(scaffold_file)
    csv_file = os.path.join(os.path.dirname(scaffold_file), os.path.basename(scaffold_file) + 'output')
    context = Context('artery csv')
    region = context.getDefaultRegion()
    region.setName('artery')

    result = region.readFile(scaffold_file)

    assert result == RESULT_OK, "Failed to load model file" + str(scaffold_file)

    # self.region = region
    field_module = region.getFieldmodule()
    with ChangeManager(field_module):
        groups = get_vessel_group_names(field_module)
        groups = get_groups_with_one_element(groups)
        lengths, radius1, radius2, relationships, distances = get_vessel_lengths(field_module, groups)
        # _write_vessel_geometry(csv_file, groups, lengths, radius1, radius2,relationships, distances)
        a=1


def get_xi_location_of_vessel(scaffold_file, csv_file):
    """

    :param scaffold_file:
    :param csv_file:
    :return:
    """
    output = os.path.join(os.path.dirname(scaffold_file), os.path.basename(scaffold_file) + 'output2')
    context = Context('artery csv')
    region = context.getDefaultRegion()

    result = region.readFile(scaffold_file)

    assert result == RESULT_OK, "Failed to load model file" + str(scaffold_file)
    field_module = region.getFieldmodule()

    with open(csv_file, 'r') as f:
        lines = f.readlines()

    new_to_old = {}
    groups = []
    data = {}
    for line in lines:
        if 'name,name' in line:
            continue
        line = line.split(',')
        new_to_old[line[1]] = line[0]
        groups.append(line[0])
        data[line[0]] = line[1:]

    parent = {}
    # for line in lines:
    #     p

    # interpolate every 1mm along the parent vessel
    cache = field_module.createFieldcache()
    coordinates = field_module.findFieldByName('coordinates').castFiniteElement()
    radius_field = field_module.findFieldByName('radius').castFiniteElement()
    # Get element iterator
    mesh = field_module.findMeshByDimension(1)

    dxi = 1.0  # mm
    elem_node_list_with_inverse_node_order = []
    node_reverse_derivative_list = []
    for group in groups:
        if group == 'Ascending aorta':
            continue
        child_name = group
        parent_name = new_to_old[data[group][5]]
        group = '"'+group+'"'  # add quotes to match the group name
        field_group = field_module.findFieldByName(group).castGroup()
        element_group = field_group.getFieldElementGroup(mesh)
        if element_group.isValid():
            mesh_group = element_group.getMeshGroup()
        else:
            print('extractPathParametersFromRegion: missing group "' + group + '"')
        elemiter = mesh_group.createElementiterator()
        element = elemiter.next()

        field_group2 = field_module.findFieldByName('"'+parent_name+'"').castGroup()
        element_group2 = field_group2.getFieldElementGroup(mesh)
        if element_group2.isValid():
            mesh_group2 = element_group2.getMeshGroup()
        else:
            print('extractPathParametersFromRegion: missing group "' + parent_name + '"')
        elemiter2 = mesh_group2.createElementiterator()
        element2 = elemiter2.next()
        while element.isValid():
            xv = []
            # for each element use Gaussian quadrature to calculate the arc length
            eft = element.getElementfieldtemplate(coordinates, 3)
            node1 = element.getNode(eft, 1)
            node2 = element.getNode(eft, 2)
            cache.setNode(node1)
            result, r1 = radius_field.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, v1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            cache.setNode(node2)
            result, r2 = radius_field.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, v2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            chv1 = v1
            chv2 = v2
            ch1id = node1.getIdentifier()
            ch2id = node2.getIdentifier()

            # get node parameters of the parent vessel
            xv = []
            # for each element use Gaussian quadrature to calculate the arc length
            eft = element2.getElementfieldtemplate(coordinates, 3)
            node1 = element2.getNode(eft, 1)
            node2 = element2.getNode(eft, 2)
            cache.setNode(node1)
            result, r1 = radius_field.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, v1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            cache.setNode(node2)
            result, r2 = radius_field.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, v2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            arcLength = getCubicHermiteArcLength(v1, d1, v2, d2)

            # point cloud
            for i in range(int(arcLength // dxi)):
                xv.append(interpolateCubicHermite(v1, d1, v2, d2, i * dxi / arcLength))

            dist1, idx1 = find_closest_point_to_vessel(chv1, xv)
            dist2, idx2 = find_closest_point_to_vessel(chv2, xv)
            if dist1 > dist2:
                # print("dist1 > dist2")
                node_reverse_derivative_list.append(ch1id)
                node_reverse_derivative_list.append(ch2id)
                elem_node_list_with_inverse_node_order.append([child_name, str(element.getIdentifier()), str(element2.getIdentifier())])
                xi = idx2 * dxi / arcLength
            else:
                xi = idx1 * dxi / arcLength
            # print(child_name, parent_name, xi, idx1, idx2, dist1, dist2)
            # import matplotlib.pyplot as plt
            # ax = plt.figure().add_subplot(projection='3d')
            # xv=np.array(xv)
            # ax.plot(xv[:, 0], xv[:, 1], xv[:, 2], 'o')
            # plt.show()
            # points[group] = xv

            element = elemiter.next()
            if element.isValid():
                print('Warning: More than one element per group')
    # for c in elem_node_list_with_inverse_node_order:
    #     print(f'{c[0]}, {c[1]}, {c[2]}')
    elem_list = [int(c[1]) for c in elem_node_list_with_inverse_node_order]
    node_reverse_derivative_list.sort()
    elem_list.sort()
    output_correct_order_of_nodes(scaffold_file, elem_list, node_reverse_derivative_list)


def output_correct_order_of_nodes(scaffold_file, elem_list, node_reverse_derivative_list):
    """
    :param scaffold_file: zinc scaffold file. arteries.exf
    :param elem_list:
    :param node_reverse_derivative_list:
    :return:
    """
    elems_started = False
    nodes_started = False
    counter = -100
    counter_node = 100
    with open(scaffold_file, 'r') as f, open(scaffold_file.replace('.exf', '_new.exf'), 'w') as g:
        for line in f:
            if line.startswith('Element template: element1'):
                elems_started = True
            elif line.startswith('Group name:'):
                elems_started = False
            if line.startswith('Node template: node1'):
                nodes_started = True
            elif line.startswith('!#mesh'):
                nodes_started = False
            if nodes_started:
                if node_reverse_derivative_list:
                    if 'Node: ' + str(node_reverse_derivative_list[0]) + '\n' in line:
                        counter_node = 0
                    if 0 < counter_node < 4:
                        line = line.strip().split()
                        g.write(f' {line[0]} {-float(line[1])}\n')
                        if counter_node == 3:
                            node_reverse_derivative_list.pop(0)
                    else:
                        g.write(line)
                    counter_node += 1
                else:
                    g.write(line)
            elif elems_started:
                if elem_list:
                    if 'Element: ' + str(elem_list[0]) + '\n' in line:
                        counter = 0
                    if counter == 2:
                        line = line.strip().split(' ')
                        g.write(f' {line[1]} {line[0]}\n')
                        elem_list.pop(0)
                    else:
                        g.write(line)
                    counter += 1
                else:
                    g.write(line)
            else:
                g.write(line)


    # get the first node of the vessel
    # get xi location of the first node of the vessel
    # find the

