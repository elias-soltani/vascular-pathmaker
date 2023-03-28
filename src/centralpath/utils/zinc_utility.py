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
    define_started = False
    with open(zinc_file, 'r') as f, open(output_zinc_mirrored, 'w') as g:
        for line in f:
            if element_started or define_started:
                g.write(line)
                continue
            if 'node2' in line:
                define_started = True
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
                dist1 = find_closest_point_to_vessel(points[group][0], points[name])
                dist2 = find_closest_point_to_vessel(points[group][-1], points[name])
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
    return dist


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
        _write_vessel_geometry(csv_file, groups, lengths, radius1, radius2,relationships, distances)
        a=1


# def my_snake_game():
#     import pygame
#     import random
#
#     # Initialize Pygame
#     pygame.init()
#
#     # Set up the game window
#     window_width = 640
#     window_height = 480
#     window = pygame.display.set_mode((window_width, window_height))
#     pygame.display.set_caption('Snake Game')
#
#     # Set up the game board
#     board_width = 20
#     board_height = 15
#     board = [[0] * board_width for i in range(board_height)]
#
#     # Set up the snake
#     snake = [(board_width // 2, board_height // 2)]
#     snake_direction = 'right'
#
#     # Set up the food
#     food = (random.randint(0, board_width - 1), random.randint(0, board_height - 1))
#
#     # Set up the game loop
#     clock = pygame.time.Clock()
#     game_over = False
#
#     while not game_over:
#         # Handle events
#         for event in pygame.event.get():
#             if event.type == pygame.QUIT:
#                 game_over = True
#             elif event.type == pygame.KEYDOWN:
#                 if event.key == pygame.K_LEFT and snake_direction != 'right':
#                     snake_direction = 'left'
#                 elif event.key == pygame.K_RIGHT and snake_direction != 'left':
#                     snake_direction = 'right'
#                 elif event.key == pygame.K_UP and snake_direction != 'down':
#                     snake_direction = 'up'
#                 elif event.key == pygame.K_DOWN and snake_direction != 'up':
#                     snake_direction = 'down'
#
#         # Move the snake
#         head_x, head_y = snake[0]
#         if snake_direction == 'left':
#             head_x -= 1
#         elif snake_direction == 'right':
#             head_x += 1
#         elif snake_direction == 'up':
#             head_y -= 1
#         elif snake_direction == 'down':
#             head_y += 1
#         snake.insert(0, (head_x, head_y))
#
#         # Check for collision with the food
#         if snake[0] == food:
#             food = (random.randint(0, board_width - 1), random.randint(0, board_height - 1))
#         else:
#             snake.pop()
#
#         # Check for collision with the walls
#         if snake[0][0] < 0 or snake[0][0] >= board_width or snake[0][1] < 0 or snake[0][1] >= board_height:
#             game_over = True
#
#         # Clear the board
#         for y in range(board_height):
#             for x in range(board_width):
#                 board[y][x] = 0
#
#         # Draw the snake
#         for x, y in snake:
#             board[y][x] = 1
#
#         # Draw the food
#         board[food[1]][food[0]] = 2
#
#         # Draw the board
#         for y in range(board_height):
#             for x in range(board_width):
#                 if board[y][x] == 0:
#                     pygame.draw.rect(window, (255, 255, 255), (x * 32, y * 32, 32, 32))
#                 elif board[y][x] == 1:
#                     pygame.draw.rect(window, (0, 255, 0), (x * 32, y * 32, 32, 32))
#                 elif board[y][x] == 2:
#                     pygame.draw.rect(window, (255, 0, 0), (x * 32, y * 32, 32, 32))
#
#         # Update the screen
#         pygame.display.update()
#
#         # Wait for the next frame
#         clock.tick(10)
#
#     # Clean up
#     pygame.quit()
#
#
# def user_snake():
#     import pygame
#     import random
#
#     # Initialize Pygame
#     pygame.init()
#
#     # Set the window size
#     window_width = 500
#     window_height = 500
#     window = pygame.display.set_mode((window_width, window_height))
#     pygame.display.set_caption("Snake Game")
#
#     # Set the colors
#     white = (255, 255, 255)
#     black = (0, 0, 0)
#     red = (255, 0, 0)
#
#     # Set the font
#     font = pygame.font.SysFont(None, 25)
#
#     # Set the clock
#     clock = pygame.time.Clock()
#
#     # Set the block size
#     block_size = 10
#
#     # Define the snake
#     def snake(block_size, snake_list):
#         for x in snake_list:
#             pygame.draw.rect(window, black, [x[0], x[1], block_size, block_size])
#
#     # Define the message
#     def message(msg, color):
#         text = font.render(msg, True, color)
#         window.blit(text, [window_width / 6, window_height / 3])
#
#     # Define the game loop
#     def gameLoop():
#         game_exit = False
#         game_over = False
#
#         # Set the initial position of the snake
#         lead_x = window_width / 2
#         lead_y = window_height / 2
#         lead_x_change = 0
#         lead_y_change = 0
#
#         # Set the initial position of the food
#         food_x = round(random.randrange(0, window_width - block_size) / 10.0) * 10.0
#         food_y = round(random.randrange(0, window_height - block_size) / 10.0) * 10.0
#
#         # Set the initial length of the snake
#         snake_list = []
#         snake_length = 1
#
#         # Game loop
#         while not game_exit:
#
#             # Game over loop
#             while game_over == True:
#                 window.fill(white)
#                 message("Game over. Press Q to quit or C to play again.", red)
#                 pygame.display.update()
#
#                 # Handle user input
#                 for event in pygame.event.get():
#                     if event.type == pygame.KEYDOWN:
#                         if event.key == pygame.K_q:
#                             game_exit = True
#                             game_over = False
#                         if event.key == pygame.K_c:
#                             gameLoop()
#
#             # Handle user input
#             for event in pygame.event.get():
#                 if event.type == pygame.QUIT:
#                     game_exit = True
#                 if event.type == pygame.KEYDOWN:
#                     if event.key == pygame.K_LEFT:
#                         lead_x_change = -block_size
#                         lead_y_change = 0
#                     elif event.key == pygame.K_RIGHT:
#                         lead_x_change = block_size
#                         lead_y_change = 0
#                     elif event.key == pygame.K_UP:
#                         lead_y_change = -block_size
#                         lead_x_change = 0
#                     elif event.key == pygame.K_DOWN:
#                         lead_y_change = block_size
#                         lead_x_change = 0
#
#             # Check if the snake hits the wall
#             if lead_x >= window_width or lead_x < 0 or lead_y >= window_height or lead_y < 0:
#                 game_over = True
#
#             # Update the position of the snake
#             lead_x += lead_x_change
#             lead_y += lead_y_change
#
#             # Draw the background
#             window.fill(white)
#
#             # Draw the food
#             pygame.draw.rect(window, red, [food_x, food_y, block_size, block_size])
#
#             # Update the snake
#             snake_head = []
#             snake_head.append(lead_x)
#             snake_head.append(lead_y)
#             snake_list.append(snake_head)
#             if len(snake_list) > snake_length:
#                 del snake_list[0]
#
#             # Check if the snake hits itself
#             for segment in snake_list[:-1]:
#                 if segment == snake_head