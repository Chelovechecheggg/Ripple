import ripple_main as ripperoni

N = 10
Vol_size = [400, 400, 360]
Vol_start = [-200, -200, -180]
coil_res = 1
vol_res = 1
plane_Xminmax = [10, 180]
plane_Zminmax = [-150, 150]
planes_Xminmax = [[10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150], [150, 170], [170, 180]]
planes_Zminmax = [[-150, -130], [-130, -110], [-110, -90], [-90, -70], [-70, -50], [-50, -30], [-30, -10], [-10, 10],
                  [10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150]]
plane_max_angle = -0.3112
plane_min_angle = 11.25
only_read = True

ripperoni.approx_N_coils('Globus2_coils', f"Globus2_coils/coils_base.txt", N, 4)
ripperoni.clone_coils('Globus2_coils', N, 16)
# ripperoni.misplace_coil(-10.3112,0,N)
# ripperoni.misplace_coil(0.3112,1,N)

coilnames = []
for i in range(int(N**2)):
    for j in range(2):
        coilnames.append(f"Globus2_coils/coil{i}_{j}.txt")
ripperoni.plot_coils(coilnames, [-15, 15], [-15, 15], [70, 100])
