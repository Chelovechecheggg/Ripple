import ripple_main as ripperoni

N = 5

#Vol_size = [400, 400, 360]
#Vol_start = [-200, -200, -180]
coil_res = 0.5
vol_res = 0.5
# plane_Xminmax = [10, 180]
# plane_Zminmax = [-150, 150]
#planes_Xminmax = [[1, 20], [20, 40], [40, 60], [60, 90]]
#planes_Zminmax = [[-100, -90], [-90, -70], [-70, -50], [-50, -30], [-30, -10], [-10, 10],
#                  [10, 30], [30, 50], [50, 70], [70, 90], [90, 100]]
planes_Xminmax = [[10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150], [150, 170], [170, 180]]
planes_Zminmax = [[-150, -130], [-130, -110], [-110, -90], [-90, -70], [-70, -50], [-50, -30], [-30, -10], [-10, 10],
                  [10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150]]
plane_max_angle = -0.3112
plane_min_angle = 11.25
only_read = True
output_filename = f'Globus3_coils/Ripple_curv_upd/Ripple_fill_5_mispl.txt'
folder_name = 'Globus3_coils'
base_coil = f"Globus3_coils/coils_base_curved_upd.txt"

ripperoni.approx_N_coils(folder_name, base_coil, N, 4)
ripperoni.clone_coils(folder_name, N, 16)
ripperoni.misplace_coil(folder_name, -0.3112, 0, N)
ripperoni.misplace_coil(folder_name, 0.3112, 1, N)


ripperoni.clear_output_file(output_filename)

for plane_Xminmax in planes_Xminmax:
    for plane_Zminmax in planes_Zminmax:
        print(plane_Xminmax, plane_Zminmax)
        ripple = ripperoni.calc_ripple(folder_name, N, vol_res, plane_Xminmax, plane_Zminmax, coil_res,
                                       plane_max_angle, plane_min_angle)
        # plot_ripple(ripple,plane_Xminmax,plane_Zminmax,vol_res)
        ripperoni.print_ripple(ripple, plane_Xminmax, plane_Zminmax, vol_res, output_filename)
