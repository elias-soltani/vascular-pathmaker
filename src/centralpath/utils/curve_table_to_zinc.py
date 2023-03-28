"""
From curve data stored in tsv file (output file of 3d Slicer), writes the elem file. The tsv file has
information about the endpoints of the curves and their lengths and radius.
"""

import csv
import os

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldFiniteElement
from opencmiss.zinc.context import Context as ZincContext
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


def read_file(input_file):
    """
    read input file
    :param input_file: tsv file. This is output of SlicerVMTK network table.
    :return: lines: List['str']
    """
    lines = []
    with open(input_file, 'r') as f:
        csv_object = csv.reader(f, delimiter='\t')
        # skip header
        next(csv_object)
        for line in csv_object:
            lines.append(line)

    return lines


def nodes_mapping(lines):
    """
    Get list of nodes from the lines. Based on the unique x,y,z.
    :param lines:
    :return:
    """
    node_map = {}
    coordinates_map = {}
    nodeIdentifier = 1
    nid_map = {}
    element_list = []
    radius = {}
    length = {}

    for line in lines:
        nodeIds = []

        x_string1 = ','.join(line[3:6])
        x_string2 = ','.join(line[6:9])
        for x_string in [x_string1, x_string2]:
            if x_string not in coordinates_map:
                nid_map[nodeIdentifier] = [float(x) for x in x_string.split(',')]
                coordinates_map[x_string] = nodeIdentifier
                nodeIdentifier += 1
                radius[coordinates_map[x_string]] = [float(line[2])]
                length[coordinates_map[x_string]] = [float(line[1])]
            else:
                radius[coordinates_map[x_string]].append(float(line[2]))
                length[coordinates_map[x_string]].append(float(line[1]))

            nodeIds.append(coordinates_map[x_string])
        element_list.append(nodeIds)

    for nid, x in nid_map.items():
        node_map[nid] = {'x': nid_map[nid], 'r': radius[nid], 'l':length[nid]}

    return node_map, element_list


def write_file(output, node_map, element_list):
    """
    Write into zinc file
    :param node_map:
    :param element_list:
    :return:
    """
    context = ZincContext("CsvToEx")
    region = context.getDefaultRegion()
    fieldmodule = region.getFieldmodule()

    coordinates = findOrCreateFieldCoordinates(fieldmodule, components_count=3)
    radius = findOrCreateFieldFiniteElement(fieldmodule, "radius", components_count=1, managed=True)
    cache = fieldmodule.createFieldcache()

    # Create nodes
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)

        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)

    mesh = fieldmodule.findMeshByDimension(1)
    # if onlyCoordinates:
    basis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
    # else:
    #     basis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eft = mesh.createElementfieldtemplate(basis)
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    for nodeIdentifier, values in node_map.items():
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, values['x'])

    elementIdentifier = 1
    for elem_nodes in element_list:
        element = mesh.createElement(elementIdentifier, elementtemplate)
        element.setNodesByIdentifier(eft, elem_nodes)
        elementIdentifier = elementIdentifier + 1

    region.writeFile(output)


def convert_tsv_to_zinc(input_file):
    """
    Get the information from tsv file (output from 3D slicer) and convert to zinc file.
    :param input_file:
    :return:
    """
    # read tsv file and parse it
    lines = read_file(input_file)
    # get list of nodes from the curves
    node_map, element_list = nodes_mapping(lines)
    # write exf
    output = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + '.exf')
    write_file(output, node_map, element_list)
