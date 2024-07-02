import ripple_main as ripperoni

N = 1.25
Vol_size = [400, 400, 360]
Vol_start = [-200, -200, -180]
coil_res = 2
vol_res = 2
plane_Xminmax = [10, 180]
plane_Zminmax = [-150, 150]
# planes_Xminmax = [[10, 30],[30, 50],[50, 70],[70, 90],[90, 110],[110, 130],[130, 150],[150, 170],[170, 180]]
# planes_Zminmax = [[-150, -130],[-130, -110],[-110, -90],[-90, -70],[-70, -50],[-50, -30],[-30,-10],[-10, 10],
#                  [10, 30],[30, 50],[50, 70],[70, 90],[90,110],[110, 130],[130, 150]]
plane_max_angle = 22.5
plane_min_angle = 11.25
output_filename = f'../Globus3_coils/Ripple.txt'
only_read = True

ripperoni.clear_output_file(output_filename)
ripperoni.approx_N_coils('../Globus3_coils', f"../Globus3_coils/coils_base_upd.txt", N, 4)
ripperoni.clone_coils('../Globus3_coils', N, 16)

ripple = ripperoni.calc_ripple('Globus3_coils', N, vol_res, plane_Xminmax, plane_Zminmax, coil_res, plane_max_angle,
                               plane_min_angle)
ripperoni.print_ripple(ripple, plane_Xminmax, plane_Zminmax, vol_res, output_filename)
