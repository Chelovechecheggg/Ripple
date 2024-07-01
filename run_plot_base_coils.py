import ripple_main as ripperoni

N = 3
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


ripperoni.approx_N_coils('Globus3_coils',f"Globus3_coils/coils_base_upd.txt", N, 4)
ripperoni.clone_coils('Globus3_coils', N, 16)
ripperoni.plot_coils([f"Globus3_coils/coil_base_1.txt",f"Globus3_coils/coil_base_2.txt",f"Globus3_coils/coil_base_3.txt",f"Globus3_coils/coil_base_4.txt",
          f"Globus3_coils/coil0_0.txt",
       f"Globus3_coils/coil1_0.txt",f"Globus3_coils/coil2_0.txt",f"Globus3_coils/coil3_0.txt",
      f"Globus3_coils/coil4_0.txt", f"Globus3_coils/coil5_0.txt", f"Globus3_coils/coil6_0.txt", f"Globus3_coils/coil7_0.txt"],
                     [10,180],[-20,20],[-150,150])
#        f"Globus3_coils/coil8_0.txt", f"Globus3_coils/coil9_0.txt", f"Globus3_coils/coil10_0.txt", f"Globus3_coils/coil11_0.txt",
#        f"Globus3_coils/coil12_0.txt", f"Globus3_coils/coil13_0.txt", f"Globus3_coils/coil14_0.txt", f"Globus3_coils/coil15_0.txt",)
