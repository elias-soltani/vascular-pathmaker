"""
utility functions to get centerlines from VMTK (Slicer) and embed to vasculatrue scaffold.
"""

from cProfile import Profile
import pstats

from utils.hermite import *
from utils.curve_table_to_zinc import *
from utils.vtk_utility import *
from utils.curve_vessel_matching import *
from utils.zinc_utility import *


if __name__ == "__main__":
    tsv_to_zic = False
    scale_radius = False
    fit_spline = False
    if tsv_to_zic:
        input_file = r'C:\Users\egha355\Desktop\work_related\Human\obj\arteries\SlicerVmtk\others\Network properties modified.tsv'
        convert_tsv_to_zinc(input_file)

    # example 2. Scale the radius by a scalar.
    if scale_radius:
        input_f = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\veins.exf'
        modify_radius_by_scalar(input_f, 15/23.1)

    # example 3. fit cubic Hermite Spline
    if fit_spline:
        vtk_input = r'C:\Users\egha355\Desktop\work_related\Human\obj\arteries\SlicerVmtk\arteries_scaledd500.vtk'
        radius_scale = 3.1
        vtk_points_to_zinc_curve(vtk_input, plot_curves=False, radius_scale=radius_scale)

        polydata = read_vtk(vtk_input, plot=False)
        curve0 = get_curve_points(polydata, 0)

        x1, x2, d1_1, d1_2 = fit_compute_derivatives(curve0, disp=True, plot=True)

    example_4 = False
    if example_4:
        input_file = r'C:\Users\egha355\Desktop\work_related\Human\obj\arteries\SlicerVmtk\veins.obj'
        vessel_file = r'C:\Users\egha355\Desktop\work_related\Human\obj\arteries\SlicerVmtk\veins.exf'
        output_name = os.path.join(os.path.dirname(input_file), 'output.exf')
        profiler = Profile()
        profiler.enable()

        groups = get_zinc_groups_from_obj_groups(vessel_file, input_file)
        write_zinc_groups(output_name, groups)

        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('tottime')
        stats.print_stats()

    example_5 = False
    if example_5:
        inputf = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\veins_out.exf'
        outputf = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\veins_out_corrected.exf'
        argon_document_output_all_zinc_groups(inputf)
        # correct_vessels_repeated_names(zinc_file)

    example_6 = False
    if example_6:
        zinc_dir = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\create group for all'
        import glob
        for file in glob.glob(zinc_dir + r'\*'):
            if os.path.isfile(file):
                mirror_zinc_file(file, [1, 0, 0, 0])

    example_7 = False
    if example_7:
        zinc_files_directory = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\create group for all'
        group_all_elements_nodes_of_zinc(zinc_files_directory)

    example_8 = False
    # outptut csv file for Finbar
    if example_8:
        output_vessel_anatomical_properties(
            r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\finbar_vessels\arteries.exf')

    example_9 = True
    if example_9:
        zfile = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\finbar_vessels\arteries_reversed_some_element_nodes_corrected.exf'
        cfile = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\finbar_vessels\arteries_v2_map_names_new.csv'
        name_map_file = r'C:\Users\egha355\Desktop\work_related\human_vasculature\cmgui\finbar_vessels\arteries_v2_map_names.csv'
        # get_xi_location_of_vessel(zfile, cfile)
        # get_children_of_vessel(cfile)
        get_radius_length_of_segments(zfile, cfile, name_map_file)

    # with open(input, 'r') as f, open(output, 'w') as g:
    #     for line in f:
    #         if 'name,length' in line:
    #             g.write(line)
    #             continue
    #         line = line.split(',')
    #         outputs = line[6].strip().split(' ')
    #         assert len(line[6].split(' ')) == int(line[7].strip())
    #         if len(outputs) > 1:
    #             cn = 1
    #             # get base name and suffix
    #             if '_L' in line[0]:
    #                 base_name = line[0].split('_L')[0]
    #                 suffix = '_L'
    #             elif '_R' in line[0]:
    #                 base_name = line[0].split('_R')[0]
    #                 suffix = '_R'
    #             else:
    #                 base_name = line[0]
    #                 suffix = ''
    #             for o in outputs:
    #                 name = ''.join([base_name, '_', str(cn), suffix])
    #                 out = ''.join([base_name, '_', str(cn + 1), suffix])
    #                 out = ' '.join([out, o])
    #                 if cn == 1:
    #                     inpt = line[5]
    #                 else:
    #                     inpt = ''.join([base_name, '_', str(cn-1), suffix])
    #                 line[0], line[5], line[6] = name, inpt, out
    #                 g.write(','.join(line))
    #                 cn += 1
    #
    #             # last output
    #             name = ''.join([base_name, '_', str(cn), suffix])
    #             inpt = ''.join([base_name, '_', str(cn-1), suffix])
    #             out = 'None'
    #             line[0], line[5], line[6] = name, inpt, out
    #             g.write(','.join(line))













    # outputs = {}
    # with open(input, 'r') as f:
    #     for line in f:
    #         if 'boundary condition type' in line:
    #             continue
    #         line = line.split(',')
    #         sp = line[5]
    #         if sp not in outputs:
    #             outputs[sp] = [line[0]]
    #         else:
    #             outputs[sp].append(line[0])
    #
    # with open(input, 'r') as f, open(output, 'w') as g:
    #     for line in f:
    #         if 'boundary condition type' in line:
    #             continue
    #         line = line.split(',')
    #         if line[0] in outputs:
    #             g.write(' '.join(outputs[line[0]])+'\n')
    #         else:
    #             g.write('None'+'\n')


    # with open(input, 'r') as f, open(output, 'w') as g:
    #     for line in f:
    #         sp = line.split(',')[6].split(' ')
    #         g.write(str(len(sp))+'\n')


