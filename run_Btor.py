import ripple_main as ripperoni

N = 3
Vol_size = [400, 400, 360]
Vol_start = [-200, -200, -180]
coil_res = 2
vol_res = 2
planes_Xminmax = [[-190, -170], [-170, -150], [-150, -130], [-130, -110], [-110, -90], [-90, -70], [-70, -50],
                  [-50, -30], [-30, -10], [-10, 10],
                  [10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150], [150, 170], [170, 190]]
planes_Yminmax = [[-190, -170], [-170, -150], [-150, -130], [-130, -110], [-110, -90], [-90, -70], [-70, -50],
                  [-50, -30], [-30, -10], [-10, 10],
                  [10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150], [150, 170], [170, 190]]
plane_max_angle = -0.3112
plane_min_angle = 11.25
only_read = True

ripperoni.approx_N_coils('Globus3_coils', f"Globus3_coils/coils_base_curved.txt", N, 4)
ripperoni.clone_coils('Globus3_coils', N, 16)
ripperoni.misplace_coil('Globus3_coils', -0.3112, 0, N)
ripperoni.misplace_coil('Globus3_coils', 0.3112, 1, N)

l = 0
for plane_Xminmax in planes_Xminmax:
    for plane_Yminmax in planes_Yminmax:
        print(plane_Xminmax, plane_Yminmax)
        ripperoni.calc_Btor(plane_Xminmax, plane_Yminmax, vol_res, coil_res, 'Globus3_coils', N, f'Btor/Btor_{l}.txt')
        l += 1
