import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import biot_savart as bs
from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.spatial.transform import Rotation as R


def parse_coil(filename):
    '''
    Parses 4 column CSV into x,y,z,I slices for coil.

    Each (x,y,z,I) entry defines a vertex on the coil.

    The current I of the vertex, defines the amount of current running through the next segment of coil, in amperes.

    i.e. (0, 0, 1, 2), (0, 1, 1, 3), (1, 1, 1, 4) means that:
    - There are 2 amps of current running between points 1 and 2
    - There are 3 amps of current running between points 2 and 3
    - The last bit of current is functionally useless.
    '''
    parsed_coil = np.loadtxt(filename, comments="#", delimiter=",", unpack=False)
    return parsed_coil.T


def plot_coils(input_filenames,xlims,ylims,zlims):
    '''
    Plots one or more coils in space.

    input_filenames: Name of the files containing the coils.
    Should be formatted appropriately.
    '''
    fig = plt.figure()
    tick_spacing = 2
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("$x$ (cm)")
    ax.set_ylabel("$y$ (cm)")
    ax.set_zlabel("$z$ (cm)")
    i = 0
    for input_filename in input_filenames:
        coil_points = np.array(parse_coil(input_filename))
        if i < 4:
            ax.plot3D(coil_points[0, :], coil_points[1, :], coil_points[2, :], lw=2, color='green')
        else:
            ax.plot3D(coil_points[0, :], coil_points[1, :], coil_points[2, :], lw=2, color='blue')
        i += 1
    ax.axes.set_xlim3d(left=xlims[0], right=xlims[1])
    ax.axes.set_ylim3d(bottom=ylims[0], top=ylims[1])
    ax.axes.set_zlim3d(bottom=zlims[0], top=zlims[1])
    plt.tight_layout()
    plt.show()


def plot_ripple(ripple_pl,Xminmax, Zminmax, vol_res):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    axes[0].set_ylabel("Z (cm)")

    for i in range(1):
        R = np.linspace(int(Xminmax[0]),int(Xminmax[1]),int((Xminmax[1] - Xminmax[0]) / vol_res)+ 1)
        Z = np.linspace(int(Zminmax[0]),int(Zminmax[1]),int((Zminmax[1] - Zminmax[0]) / vol_res)+ 1)
        contours = axes[i].contourf(R, Z, ripple_pl.T,
                                    vmin=0, vmax=ripple_pl.max(),
                                    cmap=cm.hot, levels=150)
        axes[i].set_xlabel("R (cm)")
        axes[i].set_title("Ripple")
        if i == 0:
            axes[1].set_aspect(10)
            fig.colorbar(contours, cax=axes[1], extend='both')
        plt.tight_layout()
        plt.show()


def get_plane(phi, Xmin, Xmax, Zmin, Zmax, vol_resolution):
    X = np.linspace(int(Xmin), int(Xmax), int((Xmax - Xmin) / vol_resolution) + 1)
    Z = np.linspace(int(Zmin), int(Zmax), int((Zmax - Zmin) / vol_resolution) + 1)
    plane = []
    for i in range(len(X)):
        for j in range(len(Z)):
            plane.append(np.array([X[i], 0, Z[j]]))
    plane = np.array(plane)
    r = R.from_euler('z', phi, degrees=True)
    plane = r.apply(plane)
    plane = plane.reshape((1,int((Xmax-Xmin)/vol_resolution)+1,int((Zmax-Zmin)/vol_resolution)+1,3))
    '''
    file = open('plane.txt', 'w')
    for i in plane[0]:
        for j in i:
            file.write(str(j[0])+','+str(j[1])+','+str(j[2])+'\n')
    '''

    return plane

def get_plane_Z(z, Xmin, Xmax, Ymin, Ymax, vol_resolution):
    X = np.linspace(int(Xmin), int(Xmax), int((Xmax - Xmin) / vol_resolution) + 1)
    Y = np.linspace(int(Ymin), int(Ymax), int((Ymax - Ymin) / vol_resolution) + 1)
    plane = []
    for i in range(len(X)):
        for j in range(len(Y)):
            plane.append(np.array([X[i], Y[j], z]))
    plane = np.array(plane)
    plane = plane.reshape((1,int((Xmax-Xmin)/vol_resolution)+1,int((Ymax-Ymin)/vol_resolution)+1,3))
    '''
    file = open('plane.txt', 'w')
    for i in plane[0]:
        for j in i:
            file.write(str(j[0])+','+str(j[1])+','+str(j[2])+'\n')
    '''

    return plane


def split_segments(points, N, N_base_coils):
    '''
    fig = plt.figure()
    tick_spacing = 2
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("$x$ (cm)")
    ax.set_ylabel("$y$ (cm)")
    ax.set_zlabel("$z$ (cm)")
    ax.axes.set_xlim3d(left=00, right=40)
    ax.axes.set_ylim3d(bottom=-20, top=20)
    ax.axes.set_zlim3d(bottom=120, top=160)
    '''
    splits = np.zeros(shape=(int(N_base_coils), int(N), 3))
    i = 0
    for pn in range(len(points)):
        k = 0
        segment = points[pn] - points[pn - 1]
        # print(segment)
        for sp_n in range(int(N)):
            # ax.scatter(points[pn,0], points[pn,1], points[pn,2], lw=2, color='green')
            splits[i, k] = points[pn] - segment * ((k + 1) / (int(N) + 1))
            # ax.scatter(splits[i,k,0], splits[i,k,1],splits[i,k,2], lw=2)
            k += 1
        i += 1

    # plt.tight_layout()
    # plt.show()
    return splits


def d(m, n, o, p):
    d = ((m[0] - n[0]) * (o[0] - p[0])) + ((m[1] - n[1]) * (o[1] - p[1])) + ((m[2] - n[2]) * (o[2] - p[2]))
    return d


def intersect(p1, p2, p3, p4):
    mua = (d(p1, p3, p4, p3) * d(p4, p3, p2, p1) - d(p1, p3, p2, p1) * d(p4, p3, p4, p3)) \
          / (d(p2, p1, p2, p1) * d(p4, p3, p4, p3) - d(p4, p3, p2, p1) * d(p4, p3, p2, p1))
    intersect = p1 + mua * (p2 - p1)
    return intersect


def make_approx(splits, N, N_base_coils):
    '''
    fig = plt.figure()
    tick_spacing = 2
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("$x$ (cm)")
    ax.set_ylabel("$y$ (cm)")
    ax.set_zlabel("$z$ (cm)")
    ax.axes.set_xlim3d(left=00, right=40)
    ax.axes.set_ylim3d(bottom=-20, top=20)
    ax.axes.set_zlim3d(bottom=120, top=160)
    '''

    approx = np.zeros(shape=(int(N_base_coils * (N - 1)), 3))
    j = 0
    for i in range(int(N)):
        if i == 0 or i == N - 1:
            for k in range(int(N)):
                approx[j] = intersect(splits[0, i], splits[2, -(i + 1)], splits[1, k], splits[3, -(k + 1)])

                coil_points = approx[j]
                '''
                ax.scatter(splits[0,i,0], splits[0,i,1], splits[0,i,2], lw=2,color='green')
                ax.scatter(splits[2, -(i+1), 0], splits[2, -(i+1), 1], splits[2, -(i+1), 2], lw=2, color='green')
                ax.scatter(splits[1, k, 0], splits[1, k, 1], splits[1, k, 2], lw=2, color='green')
                ax.scatter(splits[3, k, 0], splits[3, k, 1], splits[3, k, 2], lw=2, color='green')
                
                ax.scatter(coil_points[0], coil_points[1], coil_points[2], lw=2)
                '''
                j += 1
        else:
            approx[j] = intersect(splits[0, i], splits[2, -(i + 1)], splits[1, 0], splits[3, N - 1])
            '''
            coil_points = approx[j]
            ax.scatter(splits[0, i, 0], splits[0, i, 1], splits[0, i, 2], lw=2, color='green')
            ax.scatter(splits[2, i, 0], splits[2, i, 1], splits[2, i, 2], lw=2, color='green')
            ax.scatter(splits[1, N-1, 0], splits[1, N-1, 1], splits[1, N-1, 2], lw=2, color='green')
            ax.scatter(splits[3, N-1, 0], splits[3, N-1, 1], splits[3, N-1, 2], lw=2, color='green')
            ax.scatter(coil_points[0], coil_points[1], coil_points[2], lw=2)
            '''
            j += 1

            approx[j] = intersect(splits[0, i], splits[2, -(i + 1)], splits[1, N - 1], splits[3, 0])
            '''
            coil_points = approx[j]
            ax.scatter(splits[0, i, 0], splits[0, i, 1], splits[0, i, 2], lw=2, color='green')
            ax.scatter(splits[2, i, 0], splits[2, i, 1], splits[2, i, 2], lw=2, color='green')
            ax.scatter(splits[1, 0, 0], splits[1, 0, 1], splits[1, 0, 2], lw=2, color='green')
            ax.scatter(splits[3, 0, 0], splits[3, 0, 1], splits[3, 0, 2], lw=2, color='green')
            ax.scatter(coil_points[0], coil_points[1], coil_points[2], lw=2)
            '''
            j += 1

    # plt.tight_layout()
    # plt.show()

    return approx


def print_coil(path, coil, k):
    coiltxt = open(f"{path}/coil{k}.txt", 'w')
    for point in coil:
        coiltxt.write(str(point[0]) + "," + str(point[1]) + "," + str(point[2]) + "," + "1\n")
    return


def approx_N_coils(path,coil_fn, N, N_base_coils):
    all_coils = np.loadtxt(coil_fn, comments="#", delimiter="\t", unpack=False)
    # print(all_coils)
    K = int((len(all_coils) - N_base_coils + 1) / N_base_coils)
    coils = np.zeros(shape=(N_base_coils, K, 3))
    k = 0
    j = 0
    for segment in all_coils:
        if np.array_equal(segment, np.array([0, 0, 0])):
            j += 1
            k = -1
        else:
            coils[j, k] = segment
        k += 1
    # print(coils)
    points = np.zeros(shape=(N_base_coils, 3))
    approximated_coils = np.zeros(shape=(int(N_base_coils * (N - 1)), K, 3))
    for i in range(K):
        l = 0
        for ps in coils[:, i]:
            points[l] = ps
            l += 1
        splits = split_segments(points, N, N_base_coils)
        approximated_coils[:, i] = make_approx(splits, N, N_base_coils)
    print(approximated_coils)
    for k in range(len(approximated_coils)):
        print(approximated_coils[k])
        print_coil(path,approximated_coils[k], k)


def clone_coils(fpath, N, N_clones):
    for p in range(int(4 * (N - 1))):
        for l in range(N_clones):
            r = R.from_euler('z', (360/N_clones) * l, degrees=True)
            vert = np.loadtxt(f'{fpath}\coil{p}.txt', delimiter=',')
            f = open(f"{fpath}\coil{p}_{l}.txt", 'w')
            verts = []
            k = 0
            for i in vert:
                verts.append(list(np.delete(i, -1)))
            # print(verts)
            verts_rot = []
            for i in verts:
                verts_rot.append(r.apply(i))
                f.write(str(verts_rot[k][0]) + "," + str(verts_rot[k][1]) + "," + str(verts_rot[k][2]) + "," + "1\n")
                k += 1
            # print(verts_rot)
        # f.close()


def Bfield_to_Btor(Bfield,plane):
    Btor = np.zeros(shape=(len(plane),len(plane[0]),len(plane[0,0])))
    for i in range(len(plane)):
        for j in range(len(plane[i])):
            for k in range(len(plane[i,j])):
                Btor[i,j,k] = -Bfield[i, j, k, 0]*np.sin(np.arctan(plane[i,j,k,1]/plane[i,j,k,0])) + \
                       Bfield[i, j, k, 1]*np.cos(np.arctan(plane[i,j,k,1]/plane[i,j,k,0]))
    return Btor



def calc_plane_Bfield(path,N,plane,coil_res):
    #h = np.zeros(shape=(4*(N-1)*16,N_points,3))
    B_tor = np.zeros(shape=(1,len(plane[0]),len(plane[0,0])))
    B_tor2 = np.zeros(shape=(1, len(plane[0]), len(plane[0, 0])))
    for i in range(int(4*(N-1))):
        for l in range(16):

            coil = parse_coil(f"{path}/coil{i}_{l}.txt")
            sliced_coil = bs.slice_coil(coil[:3].T,coil[3:4].T,coil_res)
            B_field = bs.calculate_field(sliced_coil[0],sliced_coil[1],plane)
            B_tor2 = B_tor2 + Bfield_to_Btor(B_field,plane)
            #B_tor = B_tor + (-B_field[:, :, :, 0]*np.sin(plane_angle*np.pi/180) + B_field[:, :, :, 1]*np.cos(plane_angle*np.pi/180))
            print(i, l)
    #B_tor = B_tor.reshape(len(plane[0]),len(plane[0,0]))
    B_tor2 = B_tor2.reshape(len(plane[0]), len(plane[0, 0]))
    return B_tor2


def calc_ripple(path,N,vol_res,plane_Xminmax,plane_Zminmax,coil_res,max_angle,min_angle):
    plane_0 = get_plane(max_angle, plane_Xminmax[0], plane_Xminmax[1], plane_Zminmax[0], plane_Zminmax[1], vol_res)
    plane_1125 = get_plane(min_angle, plane_Xminmax[0], plane_Xminmax[1], plane_Zminmax[0], plane_Zminmax[1], vol_res)
    Bmax = calc_plane_Bfield(path,N,plane_0,coil_res)
    Bmin = calc_plane_Bfield(path, N, plane_1125, coil_res)
    ripple = np.zeros(shape=(len(plane_0[0]),len(plane_0[0,0])))
    ripple = (Bmax-Bmin)/(Bmax+Bmin)
    return ripple


def print_ripple(ripple, Xminmax, Zminmax, vol_res, filename):
    R = np.linspace(int(Xminmax[0]), int(Xminmax[1]), int((Xminmax[1] - Xminmax[0]) / vol_res) + 1)
    Z = np.linspace(int(Zminmax[0]), int(Zminmax[1]), int((Zminmax[1] - Zminmax[0]) / vol_res) + 1)
    file = open(filename, 'a')
    for i in range(len(R)):
        for k in range(len(Z)):
            file.write(str(R[i]) + "," + str(Z[k]) + "," + str(ripple[i,k]) +  "\n")


def misplace_coil(angle,coil_number,N):
    for i in range(int(4*(N-1))):
        coil = np.loadtxt(f'Globus3_coils/coil{i}_{coil_number}.txt', comments="#", delimiter=",", unpack=False)
        coil_coords = coil[:,:3]
        r = R.from_euler('z', angle, degrees=True)
        coil_coords = r.apply(coil_coords)
        file = open(f'Globus3_coils/coil{i}_{coil_number}.txt','w')
        for point in coil_coords:
            file.write(str(point[0]) + "," + str(point[1]) + "," + str(point[2]) + ',' + '1' "\n")
        file.close()


def calc_Btor(planeXminmax, planeYminmax, vol_res, coil_res, path, N, filename):
    planeZ = get_plane_Z(0,planeXminmax[0],planeXminmax[1],planeYminmax[0],planeYminmax[1],vol_res)
    Btor = calc_plane_Bfield(path,N,planeZ,coil_res)
    file = open(f'{path}/{filename}', 'w')
    for pointX,Bx in zip(planeZ[0],Btor):
        for pointY,By in zip(pointX,Bx):
            file.write(str(pointY[0]) + "," + str(pointY[1]) + "," + str(By) + "\n")
    file.close()


def clear_output_file(filename):
    f = open(filename, 'w')
    f.close()


# UNUSED FUNCTIONS. Archived for potential use/update in the future. None of the functions below work properly as of now
'''
def calc_volume_Bfield(N,path,only_read):
    h = []
    for i in range(4 * (N - 1)):
        for l in range(16):
            if not only_read:
                bs.write_target_volume(f"{path}/coil{i}_{l}.txt", f"{path}/coil{i}_{l}", tuple(Vol_size), tuple(Vol_start), coil_res, vol_res)

            fields, positions = bs.read_target_volume(f"{path}/coil{i}_{l}")
            # reads the volume we created
            h.append(fields)
    h_total = sum(h)
    return h_total, positions


def calc_ripple_interp(B_fields, positions, Vol_size, Vol_start, vol_res, plane_Xminmax, plane_Zminmax):
    plane_0 = get_plane(22.5, plane_Xminmax[0], plane_Xminmax[1], plane_Zminmax[0], plane_Zminmax[1], 1)
    plane_1125 = get_plane(11.25, plane_Xminmax[0], plane_Xminmax[1], plane_Zminmax[0], plane_Zminmax[1], 1)
    B_x = B_fields[:, :, :, 0]
    B_y = B_fields[:, :, :, 1]
    B_T = np.sqrt(B_x ** 2 + B_y ** 2) #This is NOT toroidal B
    x = np.linspace(Vol_start[0], Vol_size[0] + Vol_start[0], int((Vol_size[0] / vol_res) + 1))
    y = np.linspace(Vol_start[1], Vol_size[1] + Vol_start[1], int((Vol_size[1] / vol_res) + 1))
    z = np.linspace(Vol_start[2], Vol_size[2] + Vol_start[2], int((Vol_size[2] / vol_res) + 1))
    B_Ti = rgi((x, y, z), B_T)
    ripple_arr = []
    for pos_req_min, pos_req_max in zip(plane_1125, plane_0):
        pts = np.array([[pos_req_min[0], pos_req_min[1], pos_req_min[2]],
                        [pos_req_max[0], pos_req_max[1], pos_req_max[2]]])
        print(B_Ti(pts))
        if len(ripple_arr) == 48:
            print('!')
        B_min = B_Ti(pts)[0]
        B_max = B_Ti(pts)[1]
        ripple = (B_max - B_min) / (B_max + B_min)
        R = np.sqrt(pos_req_min[0] ** 2 + pos_req_min[1] ** 2)
        Z = pos_req_min[2]
        ripple_arr.append([ripple, R, Z])
    ripple_arr = np.array(ripple_arr)
    print(ripple_arr)
    return ripple_arr


def plot_Bt(Bfields, box_size, start_point, vol_resolution, which_plane='z', level=0, num_contours=50):
    
    
    Plots the set of Bfields in the given region, at the specified resolutions.

    Bfields: A 4D array of the Bfield.
    box_size: (x, y, z) dimensions of the box in cm
    start_point: (x, y, z) = (0, 0, 0) = bottom left corner position of the box AKA the offset
    vol_resolution: Division of volumetric meshgrid (generate a point every volume_resolution cm)
    which_plane: Plane to plot on, can be "x", "y" or "z"
    level : The "height" of the plane. For instance the Z = 5 plane would have a level of 5
    num_contours: THe amount of contours on the contour plot.
    
    X = np.linspace(start_point[0], box_size[0] + start_point[0], int(box_size[0] / vol_resolution) + 1)
    Y = np.linspace(start_point[1], box_size[1] + start_point[1], int(box_size[1] / vol_resolution) + 1)
    Z = np.linspace(start_point[2], box_size[2] + start_point[2], int(box_size[2] / vol_resolution) + 1)

    if which_plane == 'x':

        converted_level = np.where(X >= level)
        B_sliced = [Bfields[converted_level[0][0], :, :, i].T for i in range(3)]
        x_label, y_label = "y", "z"
        x_array, y_array = Y, Z
    elif which_plane == 'y':
        converted_level = np.where(Y >= level)
        B_sliced = [Bfields[:, converted_level[0][0], :, i].T for i in range(3)]
        x_label, y_label = "x", "z"
        x_array, y_array = X, Z
    else:
        converted_level = np.where(Z >= level)
        B_sliced = [Bfields[:, :, converted_level[0][0], i].T for i in range(3)]
        x_label, y_label = "x", "y"
        x_array, y_array = X, Y

    Bmin, Bmax = np.amin(B_sliced), np.amax(B_sliced)

    component_labels = ['T', 'y', 'z']
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    axes[0].set_ylabel(y_label + " (cm)")

    for i in range(1):
        contours = axes[i].contourf(x_array, y_array, np.sqrt(B_sliced[i] ** 2 + B_sliced[i + 1] ** 2),  #this is NOT toroidal B.
                                    vmin=Bmin, vmax=Bmax,
                                    cmap=cm.plasma, levels=num_contours)
        axes[i].set_xlabel(x_label + " (cm)")
        axes[i].set_title("$\mathcal{B}$" + "$_{}$".format(component_labels[i]))
        if i == 0:
            axes[1].set_aspect(10)
            fig.colorbar(contours, cax=axes[1], extend='both')

    plt.tight_layout()

    plt.show()


def plot_coil(*input_filenames):
    fig = plt.figure()
    tick_spacing = 2
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("$x$ (cm)")
    ax.set_ylabel("$y$ (cm)")
    ax.set_zlabel("$z$ (cm)")
    i = 0
    for input_filename in input_filenames:
        coil_points = np.array(parse_coil(input_filename))
        if i < 4:
            ax.plot3D(coil_points[0, :], coil_points[1, :], coil_points[2, :], lw=2, color='green')
        else:
            ax.plot3D(coil_points[0, :], coil_points[1, :], coil_points[2, :], lw=2, color='blue')
        i += 1
    ax.axes.set_xlim3d(left=10, right=180)
    ax.axes.set_ylim3d(bottom=-20, top=20)
    ax.axes.set_zlim3d(bottom=-150, top=150)
    plt.tight_layout()
    plt.show()
'''
