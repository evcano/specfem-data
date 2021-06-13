#!/usr/bin/env python
# This program exacts local topography from specfem3d_globe files
# Topo_File - 1 arc minute
# input:
# Topo specfem file input 
# Modeling area lat1 lon1 lat2 lon2
# Name File


import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D


def return_indexs(lat, lon, delta):
    xlat = 90.0 - lat

    if lon < 0.0:
        xlon = lon + 360.0
    if lon > 360.0:
        xlon = lon - 360.0

    return int(xlat / delta), int(xlon / delta)


def extract_topo_specfem3d_globe(filename_in, lon1, lat1, lon2, lat2, filename_out):
    # constants, do not change
    cur_dir = "/home/valeroe/packages/specfem-data/topo_bathy"
    NX_BATHY = 21600
    NY_BATHY = 10800
    dr = 1.0

    # reading file topo
    data = np.fromfile(filename_in, dtype='int16')[1:]
    print(data.shape)
    dresolution = dr / 60.0
    dx = dresolution
    dy = dresolution

    xlat1, xlon1 = return_indexs(lat1, lon1, dresolution)
    xlat2, xlon2 = return_indexs(lat2, lon2, dresolution)

    NX = abs(xlon1 - xlon2) + 1
    NY = abs(xlat1 - xlat2) + 1

    local_topo = np.zeros((NY, NX))
    xinit = xlon1
    yinit = xlat1

    if xlon1 > xlon2:
        xinit = xlon2

    if xlat1 > xlat2:
        yinit = xlat2

    yend = yinit + NY
    xend = xinit + NX

    # latitude loop
    for ilat in range(yinit, yend):
        # longitude loop
        for ilon in range(xinit, xend):
            idx = ilon + ilat * NX_BATHY
            idx2 = yend - ilat - 1
            local_topo[idx2, ilon - xinit] = float(data[idx])

    print("Size Grid (NY, NX):", local_topo.shape)
    print("lon1, lat1", lon1, lat1)
    print("lon2, lat2", lon1 + (NX-1) * dresolution, lat1 + (NY-1) * dresolution)
    print("delta:", dresolution)

    # 2d figure
    fig, ax = plt.subplots()
    ls = LightSource(azdeg=135, altdeg=45)
    rgb = ls.shade(local_topo, cmap=plt.cm.gist_earth, blend_mode='overlay')
    ima = ax.pcolormesh(local_topo, cmap=plt.cm.gist_earth)
    ax.set_aspect(aspect='equal')
    cbar = fig.colorbar(ima,
                        orientation='horizontal',
                        shrink=0.5,
                        extend='both',
                        pad=0.1,
                        format='%.f')

    ax.imshow(rgb, origin='lower')
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.set_xlabel("elevation", size=10)
    fig.savefig("%s/local_topo.png" % cur_dir, dpi=600)

    # 3d figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = np.linspace(lon1, lon1 + (NX-1) * dresolution, NX)
    y = np.linspace(lat1, lat1 + (NY-1) * dresolution, NY)
    xx, yy = np.meshgrid(x, y)
    rgb = ls.shade(local_topo,
                   cmap=cm.gist_earth,
                   vert_exag=1.0/(1000),
                   blend_mode='soft')

    ax.plot_surface(xx, yy, local_topo, facecolors=rgb)
    fig.savefig("%s/local_topo3D.png" % cur_dir, dpi=600)
    plt.show()

    # write file in specfem format
    f = open("%s/%s" % (cur_dir, filename_out), 'w')
    # latitude loop
    for j in range(0, local_topo.shape[0]):
        # longitude loop
        for i in range(0, local_topo.shape[1]):
            f.write('%8.1f\n' % (local_topo[j][i]))
    f.close()

    return NX, NY, dx, dy


if __name__ == '__main__':
    cdir = "/home/valeroe/packages/specfem-data/topo_bathy"

    filename_in = "%s/topo_bathy_etopo1_unmodified_unsmoothed.dat.bin" % cdir
    lon1 = -105.328
    lat1 = 18.0625
    lon2 = -102.033
    lat2 = 20.3217
    filename_out = "local_topo.dat"

    print("Parameters:")
    print("TopoFile:", filename_in)
    print("lon1,lat1,lon2,lat2:", lat1, lon1, lat2, lon2)

    nx, ny, dx, dy = extract_topo_specfem3d_globe(filename_in,
                                                  lon1,
                                                  lat1,
                                                  lon2,
                                                  lat2,
                                                  filename_out)
