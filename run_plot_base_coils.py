import ripple_main as ripperoni

N = 2
Vol_size = [400, 400, 360]
Vol_start = [-200, -200, -180]
coil_res = 1
vol_res = 1
plane_Xminmax = [10, 180]
plane_Zminmax = [-150, 150]
planes_Xminmax = [[10, 30],[30, 50],[50, 70],[70, 90],[90, 110],[110, 130],[130, 150],[150, 170],[170, 180]]
planes_Zminmax = [[-150, -130],[-130, -110],[-110, -90],[-90, -70],[-70, -50],[-50, -30],[-30,-10],[-10, 10],
                  [10, 30],[30, 50],[50, 70],[70, 90],[90,110],[110, 130],[130, 150]]
plane_max_angle = -0.3112
plane_min_angle = 11.25
only_read = True


ripperoni.approx_N_coils(f"Globus3_coils/coils_base_upd.txt", N, 4)
ripperoni.clone_coils('Globus3_coils', N)
ripperoni.plot_coil(f"Globus3_coils/coil_base_u_1.txt",f"Globus3_coils/coil_base_u_2.txt",f"Globus3_coils/coil_base_u_3.txt",f"Globus3_coils/coil_base_u_4.txt",
          f"Globus3_coils/coil0_0.txt",
          f"Globus3_coils/coil1_0.txt",f"Globus3_coils/coil2_0.txt",f"Globus3_coils/coil3_0.txt")
