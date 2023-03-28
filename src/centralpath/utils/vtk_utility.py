"""
utility functions for vtk binary file which is output from 3D slicer VMTK.
"""

import numpy as np
import os
import pyvista as pv

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldFiniteElement
from opencmiss.zinc.context import Context as ZincContext
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from utils.hermite import *


def read_vtk(input_file, plot=True):
    """
    Reads vtk polydata file
    :param input_file:
    :return:
    """
    # Read the VTK file using pyvista
    reader = pv.get_reader(input_file)
    # or simply polydata = pv.read(input_file)
    polydata = reader.read()

    if plot:
        polydata.plot()

    return polydata


def get_polydata_points(polydata):
    """
    Get polydata points
    :param polydata: from vtk polydata reader
    :return: points
    """
    return polydata.points


def get_radius_field(polydata):
    """
    Get polydata radius field
    :param polydata: from vtk polydata reader
    :return: radius field
    """
    return polydata.get_array('Radius')


def get_frenet_tangent_field(polydata):
    """
    Get polydata frenet tangent
    :param polydata: from vtk polydata reader
    :return: frenet tangent field
    """
    return polydata.get_array('FrenetTangent')


def get_length_field(polydata):
    """
    Get polydata length
    :param polydata: from vtk polydata reader
    :return: length field
    """
    return polydata.get_array('Length')


def get_field_by_name(polydata, name: str):
    """
    Get polydata length
    :param polydata: from vtk polydata reader
    :param name: name of the field
    :return: field
    """
    return polydata.get_array(name)


def get_curve_points(polydata, curve_index):
    """
    Get curve points
    :param polydata: from vtk polydata reader
    :param curve_index:
    :return: curve points
    """
    return polydata.cell_points(curve_index)


def get_curve_point_ids(polydata, curve_index):
    """
    Get curve point identifiers
    :param polydata: from vtk polydata reader
    :param curve_index:
    :return: curve point ids
    """
    return polydata.cell_point_ids(curve_index)


def get_number_of_curves(polydata):
    """
    Get number of curves
    :param polydata: from vtk polydata reader
    :return: number of curves
    """
    return polydata.n_cells


def get_node_element_list(polydata, method='spline', plot_curves=False, radius_scale=1.0):
    """
    extract node list and element list for using to write the linear elements.
    :param polydata: from vtk polydata reader
    :param method: either spline or line. If line it uses all the points and linear element for them. If spline, it uses
    cubic Hermite spline to compute best d1.
    :param plot_curves: If true, it plots fitted curves and data for examination purpose.
    :param radius_scale: A value to scale radius.
    :return: node list and element list
    """
    # polydata = read_vtk(input_file, plot=True)
    radius_threshold = 0.02
    # radius_scale = 1.95
    assert(method == 'spline' or method == 'line'), 'Method should be either spline or line'
    points = get_polydata_points(polydata)
    radius_field = get_radius_field(polydata)
    frenetTangent_field = get_frenet_tangent_field(polydata)
    length_field = get_length_field(polydata)
    number_of_curves = get_number_of_curves(polydata)

    node_list = []
    elem_list = []
    node_identifier = 0
    for curve_index in range(number_of_curves):
        polyline_points_ids = get_curve_point_ids(polydata, curve_index)
        curve_points = polydata.cell_points(curve_index)
        curve_points_count = len(polyline_points_ids)
        length = length_field[curve_index]

        # Clean data. 1. if the length of the curve is small skip the curve
        length_min = 5  # mm
        if length <= length_min:
            continue
        x1, x2, d1_1, d1_2 = fit_compute_derivatives(curve_points, disp=True, plot=plot_curves)

        # loop over curve points
        for i in range(curve_points_count):
            if method == 'spline':
                if 0 < i < curve_points_count - 1:
                    continue
            idx = polyline_points_ids[i]  # get its point id
            if method == 'line':
                radius = radius_field[idx]
                tangent_mag = np.linalg.norm(frenetTangent_field[idx])
                d1 = np.array([1.0, 1.0, 1.0])*length if tangent_mag < 0.00001 else \
                    frenetTangent_field[idx]*length/tangent_mag
            elif method == 'spline':
                if i == 0:
                    ir = 0
                    increment = 1
                    d1 = d1_1
                elif i == curve_points_count - 1:
                    ir = -1
                    increment = -1
                    d1 = d1_2
                # Clean data. 2. the zero radius is wrong. Get the radius from the next point.
                radius = 0
                while radius < radius_threshold:
                    radius = radius_field[polyline_points_ids[ir]]
                    ir += increment

            node_list.append({'point id': polyline_points_ids[i],
                              'x': points[idx], 'r': radius*radius_scale, 'd1': d1})
            node_identifier += 1

        if method == 'line':
            elem_list.append(polyline_points_ids)
        elif method == 'spline':
            elem_list.append([node_identifier-2, node_identifier-1])

    return node_list, elem_list


def write_endpoints_to_exnode(output, node_list, elem_list, method='spline'):
    """
    Write the endpoints extracted from each curve (line/cell) in vtk into exnode.
    :param node_list:
    :param elem_list:
    :param method:
    :return:
    """
    context = ZincContext("endpoints to exnode")
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
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.defineField(radius)
    nodetemplate.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)

    nodeIdentifier = 1
    for values in node_list:
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, values['x'].tolist())
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, values['d1'].tolist())
        radius.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, values['r'])
        nodeIdentifier += 1


    mesh = fieldmodule.findMeshByDimension(1)
    basis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
    cubicHermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

    if method == 'spline':
        eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
    else:
        eft = mesh.createElementfieldtemplate(basis)
    eftRadius = mesh.createElementfieldtemplate(basis)
    # eftRadius1.setTermNodeParameter(1, 1, 1, Node.VALUE_LABEL_VALUE, 1)
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
    result = elementtemplate.defineField(coordinates, -1, eft)
    result = elementtemplate.defineField(radius, -1, eftRadius)

    elementIdentifier = 1
    for curve in elem_list:
        for i in range(len(curve)-1):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [curve[i]+1, curve[i+1]+1])
            element.setNodesByIdentifier(eftRadius, [curve[i]+1, curve[i+1]+1])
            elementIdentifier = elementIdentifier + 1

    region.writeFile(output)


def vtk_points_to_zinc_linear(vtk_input):
    """
    Convert the vtk file to zinc file with all the points in the curves and a linear element between them.
    :param vtk_input: vtk binary file from 3D slicer VMTK.
    """
    vtk_output = os.path.join(os.path.dirname(vtk_input), os.path.basename(vtk_input).split('.')[0] + '_linear.exf')
    polydata = read_vtk(vtk_input, plot=True)
    node_list, elem_list = get_node_element_list(polydata, method='line')
    write_endpoints_to_exnode(vtk_output, node_list, elem_list)


def vtk_points_to_zinc_curve(vtk_input, plot_curves=True, radius_scale=1.0):
    """

    :param vtk_input:
    :param radius_scale: A value to scale radius.
    :return:
    """
    vtk_output = os.path.join(os.path.dirname(vtk_input), os.path.basename(vtk_input).split('.')[0] + '_curve.exf')
    polydata = read_vtk(vtk_input, plot=True)
    node_list, elem_list = get_node_element_list(polydata, method='spline', plot_curves=plot_curves,
                                                 radius_scale=radius_scale)
    write_endpoints_to_exnode(vtk_output, node_list, elem_list)


def modify_radius_by_scalar(input_f, scalar):
    """
    Scale the radius by a scalar
    :param input_f:
    :param scalar: scales the value. tested value (15/23.1)
    :return:
    """
    out_f = os.path.join(os.path.dirname(input_f), os.path.basename(input_f) + '_modified.exf')
    counter = 0
    node_started = False
    with open(input_f, 'r') as f, open(out_f, 'w') as f2:
        lines = f.readlines()
        for line in lines:
            if 'Node:' in line:
                node_started = True
                counter = 0
                f2.write(line)
            else:
                counter += 1
                if node_started and counter == 4:
                    f2.write(str(float(line)*scalar)+'\n')
                else:
                    f2.write(line)


def split_polyline(x, polydata, input_file):
    """
    Not working at the momement. It is not complete.
    :param x:
    :return:
    """
    cell_idx = 332
    polyline_points = polydata.cell_points(cell_idx)
    split_point = polydata.cell_points(333)[-1]
    # Find the indices of the points closest to the split point
    distances = np.linalg.norm(polyline_points - split_point, axis=1)
    split_idx = np.argmin(distances)

    import vtk

    # Read the VTK file
    # polydata = vtk.read('path/to/file.vtk')
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(input_file)
    reader.Update()
    # Get the polydata from the reader
    polydata = reader.GetOutput()
    # Get the points of the polydata
    points = polydata.GetPoints()
    # s=points.GetPoint(100048)

    # Get the lines of the polydata
    lines = polydata.GetLines()

    # Create a polyline splitter
    splitter = vtk.vtkSplineFilter()
    splitter.SetInputData(polydata)