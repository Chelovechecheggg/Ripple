import ripple_main as ripperoni

N = 3
Vol_size = [400, 400, 360]
Vol_start = [-200, -200, -180]
coil_res = 3
vol_res = 3
planes_Xminmax = [[10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150], [150, 170], [170, 180]]
planes_Zminmax = [[-150, -130], [-130, -110], [-110, -90], [-90, -70], [-70, -50], [-50, -30], [-30, -10], [-10, 10],
                  [10, 30], [30, 50], [50, 70], [70, 90], [90, 110], [110, 130], [130, 150]]
plane_max_angle = -0.3112
plane_min_angle = 11.25
only_read = True
output_filename = f'Globus3_coils/Ripple_test/Ripple_acc.txt'

ripperoni.approx_N_coils('Globus3_coils', f"Globus3_coils/coils_base_curved.txt", N, 4)
ripperoni.clone_coils('Globus3_coils', N, 16)
ripperoni.misplace_coil(-0.3112, 0, N)
ripperoni.misplace_coil(0.3112, 1, N)
# plot_coil(f"Globus3_coils/coil_base_u_1.txt",f"Globus3_coils/coil_base_u_2.txt",f"Globus3_coils/coil_base_u_3.txt",f"Globus3_coils/coil_base_u_4.txt",
#          f"Globus3_coils/coil0_0.txt")
#          f"Globus3_coils/coil1_0.txt",f"Globus3_coils/coil2_0.txt",f"Globus3_coils/coil3_0.txt")
#       f"Globus3_coils/coil4_0.txt", f"Globus3_coils/coil5_0.txt", f"Globus3_coils/coil6_0.txt", f"Globus3_coils/coil7_0.txt"
#        f"Globus3_coils/coil8_0.txt", f"Globus3_coils/coil9_0.txt", f"Globus3_coils/coil10_0.txt", f"Globus3_coils/coil11_0.txt",
#        f"Globus3_coils/coil12_0.txt", f"Globus3_coils/coil13_0.txt", f"Globus3_coils/coil14_0.txt", f"Globus3_coils/coil15_0.txt",)
# h_total, positions = calc_volume_Bfield(N,'Globus3_coils',only_read)
# ripple = calc_ripple_interp(h_total, positions, Vol_size, Vol_start, vol_res, plane_Xminmax, plane_Zminmax)
# plot_Bt(h_total, (60, 60, 50), (-30, -30, -25), 1, which_plane='z', level=16, num_contours=150)
# plot_ripple


# ripple = calc_ripple('Globus3_coils',N,vol_res,plane_Xminmax,plane_Zminmax,coil_res)
# print_ripple(ripple, plane_Xminmax, plane_Zminmax, vol_res,
#             f'Globus3_coils/Ripple.txt')
ripperoni.clear_output_file(output_filename)

l = 0
for plane_Xminmax in planes_Xminmax:
    for plane_Zminmax in planes_Zminmax:
        print(plane_Xminmax, plane_Zminmax)
        ripple = ripperoni.calc_ripple('Globus3_coils', N, vol_res, plane_Xminmax, plane_Zminmax, coil_res,
                                       plane_max_angle, plane_min_angle)
        # plot_ripple(ripple,plane_Xminmax,plane_Zminmax,vol_res)
        ripperoni.print_ripple(ripple, plane_Xminmax, plane_Zminmax, vol_res, output_filename)
        l += 1

'''
l = 0
for plane_Xminmax in planes_Xminmax:
    for plane_Yminmax in planes_Yminmax:
        print(plane_Xminmax,plane_Yminmax)
        ripperoni.calc_Btor(plane_Xminmax, plane_Yminmax, vol_res, coil_res, 'Globus3_coils', N,f'Btor/Btor_{l}.txt')
        l+=1


coilnames = []
for i in range(int(4*(N-1))):
    for j in range(16):
        coilnames.append(f"Globus3_coils/coil{i}_{j}.txt")
plot_coil(coilnames)
'''
