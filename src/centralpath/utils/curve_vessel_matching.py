"""
From obj and the vessel files get the name for each curve in the zinc file using the obj file. This is done by comparing
 the points on the obj groups and zinc curve points. The obj file is exported from the z anatomy .blend.
"""
import os.path

import re
import numpy as np
from scipy.spatial import cKDTree

from utils.hermite import *


def get_element_parameters(nodes, node_dict):
    """
    get element parameters
    :param nodes: a list of element node identifiers
    :param nodes: a dictionary with node identifiers as keys and {x, d1} as values
    :return: parameters
    """
    return np.array([node_dict[nodes[0]]['x'], node_dict[nodes[1]]['x'], node_dict[nodes[0]]['d1'], node_dict[nodes[1]]['d1']])


def sample_segment_points(segment_points, number_of_points):
    """
    Get random number_of_points points on the segment
    :param segment_points: a list of all points on the segment
    :param number_of_points: number of points that we want
    :return: random points on the segment
    """
    length = len(segment_points)
    if length > number_of_points:
        rng = np.random.default_rng(12345)
        random_index = rng.integers(0, length-1, size=number_of_points)
        points = []
        for index in random_index:
            points.append(segment_points[index])
        return np.array(points)
    else:
        return segment_points


def sort_points(p_s, x_c):
    """
    Sort the vessel points based on closest corresponding curve point. The points on the vessel were chosen randomly.
    :param p_s: Vessel points.
    :param x_c: Curve points.
    :return: sorted points.
    """
    closest_idx = cKDTree(p_s).query(x_c)[1]
    return p_s[closest_idx]


def get_all_groups_from_obj(input_file, translate_obj=None):
    """
    Reads all groups from the obj file and returns the group vertices.
    :param input_file: obj file, space separated, group showed by the key g and vertices by v
    :return: a dictionary of {group name, group vertices}.
    """
    translate = translate_obj
    groups = {}
    group_vertex = []
    with open(input_file, 'r') as f:
        for line in f:
            line_list = line.split()
            if line_list[0].strip() == "g":
                if group_vertex:
                    groups[group_name] = group_vertex
                group_vertex = []
                group_name = line_list[1]
            elif line_list[0].strip() == "v":
                group_vertex.append([float(line_list[i+1]) + translate[i] for i in range(3)])
        if group_vertex:
            groups[group_name] = group_vertex

    return groups


def get_nodes_parameters_with_elements_from_zinc_file(input_file):
    """
    Extract curves parameters from the zinc file. Get nodes x and d1 and elements connectivity. Assumes that nodes don't
     have versions. And elements have only two nodes.
    :param input_file: zinc file.
    :return: A dictionary of the nodes and elements with identifiers as keys and {x, d1} and node identifiers as values,
     respectively.
    """
    node_dict = {}
    elem_dict = {}
    counter = -100
    element_part = False
    with open(input_file, 'r') as f:
        reader = f.readlines()
        for line in reader:
            if 'mesh' in line:
                element_part = True

            if not element_part:
                if 'Node:' in line:
                    counter = 0
                    x = []
                    d1 = []
                    node_number = int(line.split('Node:')[1])
                else:
                    counter += 1
                    if 0 < counter < 4:
                        line_list = line.split()
                        x.append(float(line_list[0].strip()))
                        d1.append(float(line_list[1].strip()))
                    elif counter == 4:
                        node_dict[node_number] = {'x': x, 'd1': d1}
            else:
                if 'Element:' in line:
                    element_number = int(line.split('Element:')[1])
                    counter = -100
                else:
                    counter += 1
                    if counter == -98:
                        nodes_identifiers = [int(line.split()[i]) for i in range(2)]
                        elem_dict[element_number] = nodes_identifiers

    return node_dict, elem_dict


def get_closest_vessel_to_curve(vessels, curve_points):
    """
    Compare the vessels points and the curve points and find the best matching vessel name.
    :param vessels: A dictionary with vessel name as keys and vessel points as values.
    :param curve_points: A list of interpolated points on the curve.
    :return: Best matching vessel name
    """
    diff_min = 1e7
    closest_vessel_name = ''
    for vessel_name, vessel_points in vessels.items():
        sorted_points = sort_points(vessel_points, curve_points)
        diff = np.linalg.norm(sorted_points - curve_points)
        if diff <= diff_min:
            diff_min = diff
            closest_vessel_name = vessel_name

    return closest_vessel_name


def get_zinc_groups(groups, vessel_name, vessel_list_names, nodes_identifiers, element_number):
    """
    Update zinc groups with the new vessel name and the node and element groups.
    :param groups:
    :param vessel_name:
    :param vessel_list_names:
    :param nodes_identifiers:
    :param element_number:
    :return:
    """
    if vessel_name in vessel_list_names:
        if groups[vessel_name]['group_number_of_elements'] == 0:
            groups[vessel_name + f"_1"] = {
                'node group': groups[vessel_name]['node group'],
                'element group': groups[vessel_name]['element group'], 'group_number_of_elements': 0}

        groups[vessel_name]['node group'] = groups[vessel_name][
                                                     'node group'] + f', {nodes_identifiers[0]}..{nodes_identifiers[1]}'
        groups[vessel_name]['element group'] = groups[vessel_name]['element group'] + ',' + str(
            element_number)
        groups[vessel_name]['group_number_of_elements'] = groups[vessel_name]['group_number_of_elements'] + 1

        groups[vessel_name + f"_{groups[vessel_name]['group_number_of_elements'] + 1}"] = {
            'node group': f'{nodes_identifiers[0]}..{nodes_identifiers[1]}',
            'element group': str(element_number), 'group_number_of_elements': 0}
    else:
        vessel_list_names.append(vessel_name)
        groups[vessel_name] = {'node group': f'{nodes_identifiers[0]}..{nodes_identifiers[1]}',
                               'element group': str(element_number), 'group_number_of_elements': 0}

    return groups, vessel_list_names


def write_zinc_groups(output_name, groups):
    """
    Write zinc groups to the output file
    :param output_name:
    :param groups: A dictionary with group name as key and {node group, element group} as values
    :return:
    """
    # with open(os.path.join(os.path.dirname(input_f), 'output5.exf'), 'w') as f:
    with open(output_name, 'w') as f:
        for group_name, values in groups.items():
            f.write(
                f'Group name: {group_name}\n!#nodeset nodes\nNode group:\n{values["node group"]}\n'
                f'!#mesh mesh1d, dimension=1, nodeset=nodes\nElement group:\n{values["element group"]}\n')


def get_zinc_groups_from_obj_groups(zinc_file, obj_file):

    segments = get_all_groups_from_obj(obj_file, translate_obj=[0, -97, -72])
    # now we read arteries zinc file and store
    node_dict, elem_dict = get_nodes_parameters_with_elements_from_zinc_file(zinc_file)

    number_of_samples = 100
    xi = np.linspace(0, 1, number_of_samples)
    groups = {}
    segment_list_names = []
    # keep only number_of_samples points for each segment
    for segment_name, segment_points in segments.items():
        segment_points = np.array(segment_points)
        points = sample_segment_points(segment_points, number_of_samples)
        segments[segment_name] = points

    # Get the pair of the closest curve and vessel data
    for element_number, nodes_identifiers in elem_dict.items():
        parameters = get_element_parameters(nodes_identifiers, node_dict)
        x_interpolated, dx = interpolate_cubic_hermite(parameters, xi)
        min_segment_name = get_closest_vessel_to_curve(segments, x_interpolated)

        groups, segment_list_names = get_zinc_groups(groups, min_segment_name, segment_list_names, nodes_identifiers,
                                                     element_number)

    return groups


def argon_document_output_all_zinc_groups(zinc_file):
    """
    Find all the groups in zinc_file and create a argon output for the groups.
    :param zinc_file:
    :return:
    """
    dirname = os.path.dirname(zinc_file)
    basename = os.path.basename(zinc_file)
    output_argon = os.path.join(dirname, basename.split('.')[0]+'.argon')

    groups = []
    with open(zinc_file, 'r') as f:
        for line in f:
            if 'Group name:' in line:
                groups.append(line.split('Group name:')[1].strip())

    string = """
                {
                  "BoundaryMode": "ALL",
                  "CoordinateField": "coordinates",
                  "ElementFaceType": "ALL",
                  "FieldDomainType": "MESH1D",
                  "LineAttributes": {
                    "BaseSize": [
                      0,
                      0
                    ],
                    "OrientationScaleField": "radius",
                    "ScaleFactors": [
                      1,
                      1
                    ],
                    "ShapeType": "CIRCLE_EXTRUSION"
                  },
                  "Lines": {},
                  "Material": "red",
                  "RenderLineWidth": 1,
                  "RenderPointSize": 1,
                  "RenderPolygonMode": "SHADED",
                  "Scenecoordinatesystem": "LOCAL",
                  "SelectMode": "ON",
                  "SelectedMaterial": "default_selected",
                  "SubgroupField": "Right_coronary_artery_Right_coronary_artery",
                  "Tessellation": "default",
                  "Type": "LINES",
                  "VisibilityFlag": true
                },"""
    with open(output_argon, 'w') as g:
        for gro in groups:
            string1 = string.replace('"Right_coronary_artery_Right_coronary_artery"', '"\\"' + gro[1:-1] + '\\""')
            g.write(string1)


def correct_vessels_repeated_names(zinc_file):
    """
    get all the groups in the vessels and change the repeated names and make it consistent with z anatomy.
    :param zinc_file:
    :return:
    """
    dirname = os.path.dirname(zinc_file)
    basename = os.path.basename(zinc_file)
    output_zinc = os.path.join(dirname, basename.split('.')[0]+'.argon')

    check_vessel_name = True
    with open(zinc_file, 'r') as f, open(output_zinc, 'w') as g:
        for line in f:
            if 'Group name:' in line:
                line = line.split('Group name:')[1].strip()
                name = line
                if re.search(r'\.l$', line) or re.search(r'\.r$', line):
                    line = line[:-2]
                length = len(line)

                if '.l_Mesh' in line or '.r_Mesh' in line:
                    delim = '.r_' if '.r_' in line else '.l_'
                    temp = line.split(delim+'Mesh')[0]
                    name = ''.join(temp)
                    if check_vessel_name:
                        assert name + delim + 'Mesh' in line, f'{line} are not the pattern I want.{name}, {delim}'
                elif '_Mesh' in line:
                    name = line.split('_Mesh')[0]
                    if check_vessel_name:
                        assert name + '_Mesh' in line, f'{line} are not the pattern I want.{name}'
                else:
                    if '.r_' in line or '.l_' in line:
                        delim = '.r_' if '.r_' in line else '.l_'
                        name = line.split(delim)[0]
                        if check_vessel_name:
                            if not ("(M2-segment)" in line or "(M3-segment)" in line):
                                assert name + delim + name in line, f'{line} are not the pattern I want.{name}, {delim}'
                    else:
                        if re.search(r'_\d+$', line):
                            temp = line.split('_')[:-1]
                            temp = '_'.join(temp)
                            name = line[:(len(temp))//2]
                        else:
                            name = line[:length//2]
                        if check_vessel_name:
                            assert name + '_' + name in line, f'{line} are not the pattern I want.{name}'

                    if re.search(r'_\d+$', line):
                        name = name + '_' + line.split('_')[-1]
                    if '.r_' in line:
                        name = name + ' (R)'
                    elif '.l_' in line:
                        name = name + ' (L)'

                line = 'Group name: ' + '"' + name.replace('_', ' ') + '"' + '\n'

            g.write(line)