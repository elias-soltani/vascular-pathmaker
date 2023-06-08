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
from scaffoldmaker.utils.interpolation import getCubicHermiteArcLength, interpolateCubicHermite, getCubicHermiteArcLengthToXi
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier


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

    dxi = 0.1  # mm
    elem_node_list_with_inverse_node_order = []
    node_reverse_derivative_list = []
    for group in groups:
        if group == 'Ascending aorta':
            print('Ascending aorta,')
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
            xv.append(interpolateCubicHermite(v1, d1, v2, d2, 0.9999))

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
            # print(f'{child_name},{parent_name},{xi}')
            print(f'{child_name},{parent_name},{arcLength}')
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
    # elem_list = [int(c[1]) for c in elem_node_list_with_inverse_node_order]
    # node_reverse_derivative_list.sort()
    # elem_list.sort()
    # output_correct_order_of_nodes(scaffold_file, elem_list, node_reverse_derivative_list)


def get_children_of_vessel(cfile):
    """
    Gets the children of the vessels in the csv file. loops over the lines and for each parent adds the artery
     to the list of children.
    :param cfile:
    :return: list of children with their respective xi locations on the parent.
    """
    with open(cfile, 'r') as f:
        lines = f.readlines()
        children = {}
        xi_locations = {}
        arc_lengths = {}
        for line in lines:
            if line.startswith('old name,name'):
                continue
            else:
                line = line.strip().split(',')
                arc_lengths[line[1]] = float(line[2])
                if not line[6]:
                    continue
                if line[6] not in children:
                    children[line[6]] = [line[1]]
                    xi_locations[line[6]] = [line[9]]
                else:
                    children[line[6]].append(line[1])
                    xi_locations[line[6]].append(line[9])
    # sort the children based on their xi location
    for k, v in children.items():
        xi = xi_locations[k]
        children[k] = [x for _, x in sorted(zip(xi, v))]
        xi_locations[k] = sorted(xi)
    with open(cfile, 'r') as f, open(cfile.replace('.csv', '_new.csv'), 'w') as g:
        for line in f:
            line_s = line.strip().split(',')
            parent = line_s[6]
            artery = line_s[1]
            if line_s[0] == 'old name':
                continue
            if artery not in children:
                # idx_order = children[parent].index(artery)
                g.write(f'{artery}, {get_parent_name(parent, children[parent], xi_locations[parent], arc_lengths[parent], artery)}, None, {artery}, 0.0, 1.0\n')
                continue
            else:
                n = len(children[artery])
                # idx_order = children[parent].index(artery)
                if float(xi_locations[artery][0])*arc_lengths[artery] <= 5.0:
                    print(artery, 'trifurcation')
                chidx = 0
                segment = 1
                artery_name = artery
                if parent == '':
                    parent_name = 'Heart'
                else:
                    parent_name = get_parent_name(parent, children[parent], xi_locations[parent], arc_lengths[parent], artery)
                new_chidx = get_index_of_first_different_xi(xi_locations[artery], chidx+1, arc_lengths[artery])
                children_names = ' '.join([get_child_updated_name(children, xi_locations, c, arc_lengths[artery]) for c in children[artery][chidx:new_chidx]])
                if (1-float(xi_locations[artery][chidx]))*arc_lengths[artery] > 5.0:
                    children_names = get_splitted_name(artery, segment+1) + ' ' + children_names
                    artery_name = get_splitted_name(artery, segment)
                g.write(f'{artery_name}, {parent_name}, {children_names}, {artery}, 0.0, {xi_locations[artery][chidx]}\n')
                chidx = new_chidx
                while chidx < n:
                    segment += 1
                    artery_name = get_splitted_name(artery, segment)
                    parent_name = get_splitted_name(artery, segment - 1)
                    new_chidx = get_index_of_first_different_xi(xi_locations[artery], chidx+1, arc_lengths[artery])
                    children_names = ' '.join([get_child_updated_name(children, xi_locations, c, arc_lengths[artery]) for c in children[artery][chidx:new_chidx]])
                    if (1-float(xi_locations[artery][chidx]))*arc_lengths[artery] > 5.0:
                        children_names = get_splitted_name(artery, segment+1) + ' ' + children_names
                    g.write(f'{artery_name}, {parent_name}, {children_names}, {artery}, {xi_locations[artery][chidx-1]}, {xi_locations[artery][chidx]}\n')
                    chidx = new_chidx
                if (1-float(xi_locations[artery][chidx-1]))*arc_lengths[artery] > 5.0:
                    artery_name = get_splitted_name(artery, segment+1)
                    parent_name = get_splitted_name(artery, segment)
                    g.write(f'{artery_name}, {parent_name}, None, {artery}, {xi_locations[artery][chidx-1]}, 1.0\n')


def get_radius_length_of_segments(zfile, cfile, name_map_file):
    """
    Get the radius and length of the segments.
    :param zfile:
    :param cfile:
    :param name_map_file:
    :return:
    """
    context = Context('segment radius length')
    region = context.getDefaultRegion()
    result = region.readFile(zfile)

    assert result == RESULT_OK, "Failed to load model file" + str(zfile)
    field_module = region.getFieldmodule()

    with open(cfile, 'r') as f:
        data = f.readlines()

    with open(name_map_file, 'r') as f:
        lines = f.readlines()

    new_to_old = {}
    groups = []
    # data = {}
    for line in lines:
        if 'name,name' in line:
            continue
        line = line.split(',')
        new_to_old[line[1]] = line[0]
        groups.append(line[0])
        # data[line[0]] = line[1:]

    artery_whole = {}
    xi_ends = {}
    name_segments = []
    for line in data:
        line = line.split(',')
        artery_whole[line[0]] = line[3].strip()
        xi_ends[line[0]] = [float(line[4]), float(line[5])]
        name_segments.append(line[0])

    cache = field_module.createFieldcache()
    coordinates = field_module.findFieldByName('coordinates').castFiniteElement()
    radius_field = field_module.findFieldByName('radius').castFiniteElement()
    # Get element iterator
    mesh = field_module.findMeshByDimension(1)

    for segment in name_segments:
        artery = artery_whole[segment]
        xi_start = xi_ends[segment][0]
        xi_end = xi_ends[segment][1]
        group = '"' + new_to_old[artery] + '"'  # add quotes to match the group name
        field_group = field_module.findFieldByName(group).castGroup()
        element_group = field_group.getFieldElementGroup(mesh)
        if element_group.isValid():
            mesh_group = element_group.getMeshGroup()
        else:
            print('extractPathParametersFromRegion: missing group "' + group + '"')
        elemiter = mesh_group.createElementiterator()
        element = elemiter.next()
        while element.isValid():
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
            length = getCubicHermiteArcLengthToXi(v1, d1, v2, d2, xi_end) - getCubicHermiteArcLengthToXi(v1, d1, v2, d2, xi_start)
            radius1 = (1-xi_start)*r1 + xi_start*r2
            radius2 = (1-xi_end)*r1 + xi_end*r2
            print(f'{segment},{length:.3f},{radius1:.3f},{radius2:.3f}')
            element = elemiter.next()
            if element.isValid():
                print('Warning: More than one element per group')


def get_parent_name(parent, parent_children, parent_xi, parent_length, artery):
    """
    Get the segmented parent name of the artery. If the parent is not segmented, then return the parent name.
    :param parent: parent original name
    :param parent_children:
    :param parent_xi:
    :param parent_length:
    :param artery: artery name
    :return:
    """
    if (1-float(parent_xi[0]))*parent_length < 5.0:
        return parent
    else:
        return get_splitted_name(parent, get_vessel_segment_number_on_parent(parent_children, artery, parent_xi, parent_length))


def get_index_of_first_different_xi(xi, idx, arc_length):
    """
    Get the index of the first xi value that is different from the previous one. i.e. xi[j] - xi[i] >= 0.05
    :param xi: list of xi values.
    :param idx: index of current artery segment.
    :param arc_length: arc length of the current artery segment.
    :return:
    """
    j = idx
    while j < len(xi) and (float(xi[j]) - float(xi[idx-1]))*arc_length <= 5.0:
        j += 1
    return j


def get_vessel_segment_number_on_parent(parent_children, name,  parent_xi, parent_length):
    """
    Get the vessel segment number on the parent.
    :param parent_children:
    :param name:
    :param parent_xi:
    :param parent_length:
    :return:
    """
    idx_order = parent_children.index(name)
    if idx_order == 0:
        return 1
    i, j = 1, 1
    while i <= idx_order:
        if (float(parent_xi[i]) - float(parent_xi[i-1]))*parent_length >= 5.0:
            j += 1
        i += 1
    return j


def get_splitted_name(name, i):
    """
    :param name: name of the vessel
    :param i: index of the vessel
    :return: name of the vessel with the index appended to the end
    """
    if name.endswith('_L'):
        return name[:-2] + '_' + str(i) + '_L'
    elif name.endswith('_R'):
        return name[:-2] + '_' + str(i) + '_R'
    else:
        return name + '_' + str(i)


def get_child_updated_name(children, xi_locations, name, arc_length):
    """
    return child name with the correct index. If the child has at least one child and the child is not at xi>=0.97 then
    returns get_splitted_name(child, 1)
    :param children:
    :param xi_locations:
    :param name:
    :param i:
    :return:
    """
    if name not in children:
        return name
    n = len(children[name])
    if n == 0:
        return name
    # elif n > 1:
    #     return get_splitted_name(name, 1)
    else:
        if (1-float(xi_locations[name][0]))*arc_length <= 5.0:
            return name
        else:
            return get_splitted_name(name, 1)


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


def get_terminals(cfile):
    """
    Get the terminals of the vessels. Read the cfile and if the vessel output is None, then make a terminal.
    :param cfile:
    :return:
    """
    with open(cfile, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('artery segment'):
            continue
        line = line.strip().split(',')
        children = line[2].strip().split(' ')
        if children[0] == 'None':
            print(f'{line[0]}_T,,,,,,,,,terminal, pp')
        # else:
        # if len(children) > 1:
        #     print('split_junction, pv')
        # else:
        #     print('arterial, pv')


def get_maximum_number(zincfile):
    """
    Get the maximum number of nodes, lines, faces and elements identifiers of the zinc file.
    :param zincfile:
    :return: maximum node, line, face and element identifiers
    """
    context = Context("maximum_number")
    region = context.getDefaultRegion()
    region.readFile(zincfile)
    fieldmodule = region.getFieldmodule()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    mesh1d = fieldmodule.findMeshByDimension(1)
    mesh2d = fieldmodule.findMeshByDimension(2)
    mesh3d = fieldmodule.findMeshByDimension(3)
    nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes))
    lineIdentifier = max(1, getMaximumElementIdentifier(mesh1d))
    faceIdentifier = max(1, getMaximumElementIdentifier(mesh2d))
    elementIdentifier = max(1, getMaximumElementIdentifier(mesh3d))

    return [nodeIdentifier, lineIdentifier, faceIdentifier, elementIdentifier]


def get_coordinates_name(zincfile):
    """
    Get the coordinates field name of the zinc file.
    :param zincfile:
    :return: coordinates field name
    """
    context = Context("coordinates_name")
    region = context.getDefaultRegion()
    region.readFile(zincfile)
    fieldmodule = region.getFieldmodule()
    if fieldmodule.findFieldByName('coordinates').isValid():
        return 'coordinates'
    elif fieldmodule.findFieldByName('fitted coordinates').isValid():
        return 'fitted coordinates'


def get_zinc_files_in_directory(directory):
    """
    Get the zinc files in the directory.
    :param directory:
    :return: list of zinc files in the directory
    """
    zinc_files = []
    for file in os.listdir(directory):
        if file.endswith(".exf"):
            zinc_files.append(os.path.join(directory, file))
    return zinc_files


def find_files_with_coordinates_name_coordinates(directory):
    """
    Find the files that have coordinates field name as coordinates.
    :param directory:
    :return: list of zinc files that have coordinates field name as coordinates
    """
    zinc_files = get_zinc_files_in_directory(directory)
    zinc_files_with_coordinates_name_coordinates = []
    for zinc_file in zinc_files:
        if get_coordinates_name(zinc_file) == 'coordinates':
            zinc_files_with_coordinates_name_coordinates.append(os.path.basename(zinc_file))
    return zinc_files_with_coordinates_name_coordinates


def combine_files(directory):
    """
    Combine the files in the directory into one file.
    :param directory:
    :param output_file:
    :return:
    """
    context = Context("combine_files")
    region = context.getDefaultRegion()
    zinc_files = get_zinc_files_in_directory(directory)
    for zinc_file in zinc_files:
        region.readFile(zinc_file)
    region.writeFile(os.path.join(directory, 'zinc.combined'))


def generate_cmgui_commands_to_write_combined_file(directory):
    """
    Generate cmgui commands to write the combined file.
    :param directory:
    :return: cmgui command file to write the combined file
    """
    zinc_files = get_zinc_files_in_directory(directory)
    output = os.path.join(directory, 'zinc_commands.combined')
    identifiers = np.array([0, 0, 0, 0])
    with open(output, 'w') as g:
        for zinc_file in zinc_files:
            g.write(f'gfx read elements {zinc_file} node_offset {identifiers[0]} line_offset {identifiers[1]}'+
                    f' face_offset {identifiers[2]} element_offset {identifiers[3]}' + '\n')
            identifiers += np.array(get_maximum_number(zinc_file))
        g.write('gfx write elem node combined.ex' '\n')


{
    "BoundaryMode": "BOUNDARY",
    "CoordinateField": "fitted coordinates",
    "ElementFaceType": "ALL",
    "FieldDomainType": "MESH2D",
    "Material": "bone",
    "RenderLineWidth": 1,
    "RenderPointSize": 1,
    "RenderPolygonMode": "SHADED",
    "Scenecoordinatesystem": "LOCAL",
    "SelectMode": "ON",
    "SelectedMaterial": "default_selected",
    "SubgroupField": "tibia_left",
    "Surfaces": {},
    "Tessellation": "default",
    "Type": "SURFACES",
    "VisibilityFlag": true
}