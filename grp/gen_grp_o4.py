#!/usr/bin/env python3

import os
from mpmath import mp
import numpy as np
from quat import Quat
from grp_o4 import GrpElemO4

mp.dps = 30

sqrt1_2 = 1 / mp.sqrt(2)
alpha = mp.cos(mp.pi / 5)
beta = mp.cos(2 * mp.pi / 5)


def FixElem(g):
    # apply sign convention to avoid double counting degenerate quaternions
    lw = mp.chop(g.ql.w)
    lx = mp.chop(g.ql.x)
    ly = mp.chop(g.ql.y)
    lz = mp.chop(g.ql.z)
    rw = mp.chop(g.qr.w)
    rx = mp.chop(g.qr.x)
    ry = mp.chop(g.qr.y)
    rz = mp.chop(g.qr.z)

    s = 1
    if lw < 0:
        s = -1
    elif lw == 0 and lx < 0:
        s = -1
    elif lw == 0 and lx == 0 and ly < 0:
        s = -1
    elif lw == 0 and lx == 0 and ly == 0 and lz < 0:
        s = -1

    ql = Quat(lw * s, lx * s, ly * s, lz * s)
    qr = Quat(rw * s, rx * s, ry * s, rz * s)
    return GrpElemO4(ql, qr, g.inv)


def ElemName(g):
    # make a string to identify this quaternion
    lw = mp.nstr(mp.chop(g.ql.w), 20)
    lx = mp.nstr(mp.chop(g.ql.x), 20)
    ly = mp.nstr(mp.chop(g.ql.y), 20)
    lz = mp.nstr(mp.chop(g.ql.z), 20)
    rw = mp.nstr(mp.chop(g.qr.w), 20)
    rx = mp.nstr(mp.chop(g.qr.x), 20)
    ry = mp.nstr(mp.chop(g.qr.y), 20)
    rz = mp.nstr(mp.chop(g.qr.z), 20)

    return ','.join([str(g.inv), lw, lx, ly, lz, rw, rx, ry, rz])


def VertexName(v):
    # make a string to identify this vertex
    x = mp.nstr(mp.chop(v[0]), 20)
    y = mp.nstr(mp.chop(v[1]), 20)
    z = mp.nstr(mp.chop(v[2]), 20)
    w = mp.nstr(mp.chop(v[3]), 20)

    return ','.join([x, y, z, w])


def GenerateGroup(q):

    print(f'Generating group o4q{q}...')

    # make a list of group generators
    gen = []
    if q == 3:
        gen.append(GrpElemO4(
            Quat(0, 0, sqrt1_2, sqrt1_2),
            Quat(0, 0, -sqrt1_2, -sqrt1_2), -1))
        gen.append(GrpElemO4(
            Quat(beta, -0.5, -alpha, 0),
            Quat(-alpha, 0.5, -beta, 0), 1))
    elif q == 4:
        gen.append(GrpElemO4(
            Quat(sqrt1_2, sqrt1_2, 0, 0),
            Quat(sqrt1_2, sqrt1_2, 0, 0), 1))
        gen.append(GrpElemO4(
            Quat(0.5, 0.5, 0.5, 0.5),
            Quat(-0.5, 0.5, 0.5, 0.5), -1))
    elif q == 5:
        gen.append(GrpElemO4(
            Quat(0.5, alpha, beta, 0),
            Quat(0.5, alpha, beta, 0), 1))
        gen.append(GrpElemO4(
            Quat(alpha, beta, 0.5, 0.0),
            Quat(1.0, 0.0, 0.0, 0.0), -1))

    # start with the identity element only
    G = [GrpElemO4.Identity()]
    hash_table = {ElemName(G[0]): 1}

    while (True):
        # new elements
        new_ones = []

        # multiply all current elements times all generators
        for g in G:
            for h in gen:
                # generate a new group element
                gh = FixElem(g * h)

                # skip if it's not new
                hash_string = ElemName(gh)
                if hash_string in hash_table:
                    continue

                # add it to the list of new ones
                new_ones.append(gh)
                hash_table[hash_string] = 1

        # exit when no new elements were generated
        if (len(new_ones) == 0):
            break

        # add new elements to the list
        G += new_ones

    print(f'Group has {len(G)} elements')

    # generate a set of vertices
    vertices = []
    vertex_table = {}

    z_hat = np.array([0.0, 0.0, 0.0, 1.0])

    for g in G:
        v = g * z_hat
        hash_string = VertexName(v)
        if hash_string in vertex_table:
            continue

        vertices.append(v)
        vertex_table[hash_string] = 1

    # determine the edge length
    edge_length = 2.0
    for s in range(1, len(vertices)):
        length = np.linalg.norm(vertices[s] - vertices[0])
        if length < edge_length:
            edge_length = length

    # generate a set of faces
    cells = []
    for s1 in range(len(vertices)):
        for s2 in range(s1 + 1, len(vertices)):
            length12 = np.linalg.norm(vertices[s1] - vertices[s2])
            if not mp.almosteq(length12, edge_length):
                continue
            for s3 in range(s2 + 1, len(vertices)):
                length13 = np.linalg.norm(vertices[s1] - vertices[s3])
                if not mp.almosteq(length13, edge_length):
                    continue
                length23 = np.linalg.norm(vertices[s2] - vertices[s3])
                if not mp.almosteq(length23, edge_length):
                    continue

                for s4 in range(s3 + 1, len(vertices)):
                    length14 = np.linalg.norm(vertices[s1] - vertices[s4])
                    if not mp.almosteq(length14, edge_length):
                        continue
                    length24 = np.linalg.norm(vertices[s2] - vertices[s4])
                    if not mp.almosteq(length24, edge_length):
                        continue
                    length34 = np.linalg.norm(vertices[s3] - vertices[s4])
                    if not mp.almosteq(length34, edge_length):
                        continue

                    cells.append([s1, s2, s3, s4])

    print(f'Polytope has {len(vertices)} vertices and {len(cells)} cells')
    print('Writing group data to file...\n')

    # write group elements to a data file
    file = open(f'elem/o4q{q}.dat', 'w')
    for i, g in enumerate(G):
        lw = mp.chop(g.ql.w)
        lx = mp.chop(g.ql.x)
        ly = mp.chop(g.ql.y)
        lz = mp.chop(g.ql.z)
        rw = mp.chop(g.qr.w)
        rx = mp.chop(g.qr.x)
        ry = mp.chop(g.qr.y)
        rz = mp.chop(g.qr.z)
        file.write('%d %+d' % (i, g.inv))
        file.write(' %+.18f %+.18f %+.18f %+.18f' % (lw, lx, ly, lz))
        file.write(' %+.18f %+.18f %+.18f %+.18f\n' % (rw, rx, ry, rz))
    file.close()

    # write polytope lattice to a data file
    file = open(f'lattice/o4q{q}.dat', 'w')
    file.write('begin_sites\n')
    file.write(f'n_sites {len(vertices)}\n')
    for i, v in enumerate(vertices):
        theta = 0.0
        phi = 0.0
        cos_psi = mp.chop(v[3])
        if mp.almosteq(cos_psi, 1.0):
            psi = 0.0
        elif mp.almosteq(cos_psi, -1.0):
            psi = mp.pi
        else:
            psi = mp.acos(cos_psi)
            cos_theta = mp.chop(v[2]) / mp.sin(psi)
            if mp.almosteq(cos_theta, 1.0):
                theta = 0.0
            elif mp.almosteq(cos_theta, -1.0):
                theta = mp.pi
            else:
                theta = mp.acos(cos_theta)
            phi = mp.atan2(mp.chop(v[1]), mp.chop(v[0]))

        psi = mp.chop(psi)
        theta = mp.chop(theta)
        phi = mp.chop(phi)
        file.write('%d 1.0 0 %+.18f %+.18f %+.18f\n' % (i, psi, theta, phi))
    file.write('end_sites\n')

    file.write('begin_cells\n')
    file.write(f'n_cells {len(cells)}\n')
    for i, c in enumerate(cells):
        file.write('%d 1.0 %d %d %d %d\n' % (i, c[0], c[1], c[2], c[3]))
    file.write('end_cells\n')
    file.close()


if __name__ == '__main__':
    os.makedirs('elem', exist_ok=True)
    os.makedirs('lattice', exist_ok=True)
    os.makedirs('site_g', exist_ok=True)
    GenerateGroup(3)
    GenerateGroup(4)
    GenerateGroup(5)
