import plotly.graph_objects as go
import meshplot
import numpy as np
from pyscf import lib
from pyscf.dft import numint, gen_grid
from skimage import measure
import meshplot as mp
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d


import imageio
from skimage import measure
from pyscf import gto, scf
from pyscf.tools import cubegen


def minmaxnorm(x, minx, maxx):
    return np.divide(x - minx, maxx - minx)


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


def density(mol, outfile, dm, nx=80, ny=80, nz=80):
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
    return coords



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
           [ 'H'  , ( 2.497289 ,   6.120281 , 0.00  ) ]]



naphthalene = [['C', (2.423247, 0.705674, 0.000000)],
            ['C', (1.239335, 1.396987, 0.000000  )],
            ['C', (0.000000, 0.713185, 0.000000  )],
            ['C', (0.000000, -0.713185, 0.000000 )],
            ['C', (1.239335, -1.396987, 0.000000 )],
            ['C', (2.423247, -0.705674, 0.000000 )],
            ['H', (-1.237376, 2.476742, 0.000000 )],
            ['H', (3.361100, 1.238312, 0.000000  )],
            ['H', (1.237376, 2.476742, 0.000000  )],
            ['C', (-1.239335, 1.396987, 0.000000 )],
            ['C', (-1.239335, -1.396987, 0.000000)],
            ['H', (1.237376, -2.476742, 0.000000 )],
            ['H', (3.361100, -1.238312, 0.000000 )],
            ['C', (-2.423247, -0.705674, 0.000000)],
            ['C', (-2.423247, 0.705674, 0.000000 )],
            ['H', (-1.237376, -2.476742, 0.000000)],
            ['H', (-3.361100, -1.238312, 0.000000)],
            ['H', (-3.361100, 1.238312, 0.000000 )]]


atoms_x = []
atoms_y = []
atoms_z = []
for atom in naphthalene:
    atoms_x.append(atom[1][0])
    atoms_y.append(atom[1][1])
    atoms_z.append(atom[1][2])
print(min(atoms_x), min(atoms_y))
print(max(atoms_x), max(atoms_y))

atoms_x_n = []
atoms_y_n = []
atoms_z_n = []

for x, y, z in zip(atoms_x, atoms_y, atoms_z):
    atoms_x_n.append(minmaxnorm(x, min(atoms_x), max(atoms_x)))
    atoms_y_n.append(minmaxnorm(y, min(atoms_y), max(atoms_y)))
    atoms_z_n.append(minmaxnorm(z, min(atoms_z), max(atoms_z)))



mol = gto.M(atom=naphthalene)
mf = scf.RHF(mol)
mf.scf()
cord = cubegen.density(mol, "h2.cube", mf.make_rdm1())


cube = parse_cube('h2.cube')



contours = []
for i in range(0, 80):
    if i % 5 == 0:
        fig = plt.figure(figsize=(10, 10))
        plt.contourf(cube["data"][:, :, i])
        plt.savefig(f"img/contour{i}.png")
        contours.append(f"img/contour{i}.png")

conts = []
for filename in contours:
    conts.append(imageio.imread(filename))
imageio.mimsave('img/contour_animation.gif', conts)

pot_range = np.linspace(0.08, 0.3, 20)
i = 0
files = []
camera = dict(
    eye=dict(x=0., y=0, z=1.5)
)
for coff in pot_range:
    vert, faces, norm, values = measure.marching_cubes(cube['data'], coff, spacing=(0.1, 0.1, 0.1))

    fig = go.Figure(data=[
        go.Mesh3d(
            x=vert[:, 0],
            y=vert[:, 1],
            z=vert[:, 2],
            colorbar_title='z',
            colorscale=[[0, 'gold'],
                        [0.5, 'mediumturquoise'],
                        [1, 'magenta']],
            # i, j and k give the vertices of triangles
            intensity=[1.0, 2.0],
            i=faces[:, 0],
            j=faces[:, 1],
            k=faces[:, 2],
            name='y',
            showscale=True,
            opacity=0.8
        )
    ])
    fig.update_layout(
        scene=dict(
            xaxis=dict(nticks=4, range=[0, 7.5], ),
            yaxis=dict(nticks=4, range=[0, 7.5], ),
            zaxis=dict(nticks=4, range=[0, 7.5], ), ),
        width=700,
        scene_camera=camera,
    )

    fig.update_scenes(xaxis_visible=False, yaxis_visible=False,zaxis_visible=False )

    fig.write_image(f"img/{i}.png")
    files.append(f"img/{i}.png")
    i = i + 1
    print(i)

images = []
for filename in files:
    images.append(imageio.imread(filename))
imageio.mimsave('img/animation.gif', images)