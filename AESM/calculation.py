#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#
import meshplot
import numpy as np
from pyscf import lib
from pyscf.dft import numint, gen_grid
from skimage import measure
import meshplot as mp
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d
import webbrowser
from html2image import Html2Image
import imageio
meshplot.offline()

"""
Gaussian cube file format
"""


def density(mol, outfile, dm, nx=100, ny=100, nz=100):
    coord = mol.atom_coords()
    box = np.max(coord, axis=0) - np.min(coord, axis=0) + 4
    boxorig = np.min(coord, axis=0) - 2
    xs = np.arange(nx) * (box[0] / nx)
    ys = np.arange(ny) * (box[1] / ny)
    zs = np.arange(nz) * (box[2] / nz)
    coords = lib.cartesian_prod([xs, ys, zs])
    coords = np.asarray(coords, order="C") - (-boxorig)

    nao = mol.nao_nr()
    ngrids = nx * ny * nz
    blksize = min(200, ngrids)
    rho = np.empty(ngrids)
    for ip0, ip1 in gen_grid.prange(0, ngrids, blksize):
        ao = numint.eval_ao(mol, coords[ip0:ip1])
        rho[ip0:ip1] = numint.eval_rho(mol, ao, dm)
    rho = rho.reshape(nx, ny, nz)

    with open(outfile, "w") as f:
        f.write("Density in real space\n")
        f.write("Comment line\n")
        f.write("]" % mol.natm)
        f.write(" .8f .8f .8f\n" % tuple(boxorig.tolist()))
        f.write("] .8f .8f .8f\n" % (nx, xs[1], 0, 0))
        f.write("] .8f .8f .8f\n" % (ny, 0, ys[1], 0))
        f.write("] .8f .8f .8f\n" % (nz, 0, 0, zs[1]))
        for ia in range(mol.natm):
            chg = mol.atom_charge(ia)
            f.write("%5d %f" % (chg, chg))
            f.write(" .8f .8f .8f\n" % tuple(coord[ia]))
        fmt = " .8e" * nz + "\n"
        for ix in range(nx):
            for iy in range(ny):
                f.write(fmt % tuple(rho[ix, iy].tolist()))

def parse_cube(filename):
    """ Parses a cube file, returning a dict of the information contained.
        The cubefile itself is stored in a np array. """
    with open(filename) as fp:
        results = {}

        # skip over the title
        fp.readline()
        fp.readline()

        origin = fp.readline().split()
        natoms = int(origin[0])
        results['minx'] = minx = float(origin[1])
        results['miny'] = miny = float(origin[2])
        results['minz'] = minz = float(origin[3])

        infox = fp.readline().split()
        numx = int(infox[0])
        incx = float(infox[1])
        results['incx'] = incx
        results['numx'] = numx
        results['maxx'] = minx + incx * numx

        infoy = fp.readline().split()
        numy = int(infoy[0])
        incy = float(infoy[2])
        results['incy'] = incy
        results['numy'] = numy
        results['maxy'] = miny + incy * numy

        infoz = fp.readline().split()
        numz = int(infoz[0])
        incz = float(infoz[3])
        results['incz'] = incz
        results['numz'] = numz
        results['maxz'] = minz + incz * numz

        atnums = []
        coords = []
        for atom in range(natoms):
            coordinfo = fp.readline().split()
            atnums.append(int(coordinfo[0]))
            coords.append(list(map(float, coordinfo[2:])))
        results['atom_numbers'] = np.array(atnums)
        results['atom_coords'] = np.array(coords)

        data = np.array([ float(entry) for line in fp for entry in line.split() ])
        if len(data) != numx*numy*numz:
            raise Exception("Amount of parsed data is inconsistent with header in Cube file!")
        results['data'] = data.reshape((numx,numy,numz))

        return results

if __name__ == "__main__":
    from pyscf import gto, scf
    from pyscf.tools import cubegen
    from selenium import webdriver
    from webdriver_manager.firefox import GeckoDriverManager



    shading = {"flat": True,  # Flat or smooth shading of triangles
               "wireframe": True, "wire_width": 0.03, "wire_color": "black",  # Wireframe rendering
               "width": 1000, "height": 1000,  # Size of the viewer canvas
               "antialias": True,  # Antialising, might not work on all GPUs
               "scale": 2.0,  # Scaling of the model
               "side": "DoubleSide",  # FrontSide, BackSide or DoubleSide rendering of the triangles
               "colormap": None, None: [None, None],  # Colormap and normalization for colors
               "background": "#ffffff",  # Background color of the canvas
               "line_width": 1.0, "line_color": "black",  # Line properties of overlay lines
               "bbox": False,  # Enable plotting of bounding box
               "point_color": "white", "point_size": 0.15  # Point properties of overlay points
               }
    benzene = [[ 'C'  , ( 4.673795 ,   6.280948 , 0.00  ) ],
           [ 'C'  , ( 5.901190 ,   5.572311 , 0.00  ) ],
           [ 'C'  , ( 5.901190 ,   4.155037 , 0.00  ) ],
           [ 'C'  , ( 4.673795 ,   3.446400 , 0.00  ) ],
           [ 'C'  , ( 3.446400 ,   4.155037 , 0.00  ) ],
           [ 'C'  , ( 3.446400 ,   5.572311 , 0.00  ) ],
           [ 'H'  , ( 4.673795 ,   7.376888 , 0.00  ) ],
           [ 'H'  , ( 6.850301 ,   6.120281 , 0.00  ) ],
           [ 'H'  , ( 6.850301 ,   3.607068 , 0.00  ) ],
           [ 'H'  , ( 4.673795 ,   2.350461 , 0.00  ) ],
           [ 'H'  , ( 2.497289 ,   3.607068 , 0.00  ) ],
           [ 'N'  , ( 2.497289 ,   6.120281 , 0.00  ) ]]

    mol = gto.M(atom=benzene)
    mf = scf.RHF(mol)
    mf.scf()
    cubegen.density(mol, "h2.cube", mf.make_rdm1())
    cube = parse_cube('h2.cube')
    print(cube['data'].min(), cube['data'].max())

    rng = np.linspace(0.03, 0.35, 10)
    i = 0
    files = []
    for cutoff in rng:
        vert, faces, norm, values = measure.marching_cubes(cube['data'],
                                                           cutoff,
                                                           spacing=(0.1, 0.1, 0.1))
        vert = vert
        print(np.shape(vert))
        exit()
        p = mp.plot(vert, faces, c=norm[:, 0],
                    shading=shading, return_plot=True,
                    filename=f'test{np.round(cutoff, 2)}.html')
        driver = webdriver.Firefox(executable_path=GeckoDriverManager().install())

        driver.get(f'file:///home/lr0n/Documents/Uni/AESM/test{np.round(cutoff, 2)}.html')
        time.sleep(1)
        driver.get_screenshot_as_file(f"img/screenshot{i}.png")
        files.append(f"img/screenshot{i}.png")
        driver.quit()
        i = i + 1
        #webbrowser.open(f'test{np.round(cutoff, 2)}.html', new=2)  # open in new tab

    print(files)
    images = []
    for filename in files:
        images.append(imageio.imread(filename))
    imageio.mimsave('img/animation.gif', images)







