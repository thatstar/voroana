from itertools import product
import numpy as np
from scipy.special import sph_harm
from sympy.physics.wigner import wigner_3j


__all__ = ['BondOrientOrder']


def _wigner_3j_values(l):
    m = range(-l, l + 1)
    inds = filter(lambda x: sum(x) == 0, product(m, m, m))
    res = {i: float(wigner_3j(l, l, l, *i).evalf(20)) for i in inds}

    return res


class BondOrientOrder:
    def __init__(self, atoms, voro_info, weighted=False):
        self.n = len(atoms)
        assert self.n == len(voro_info['neighbors'])
        self.positions = atoms.positions
        self.neighbors = voro_info['neighbors']
        self.neighbor_vectors = voro_info['neighbor_vectors']
        self.weighted = weighted
        if self.weighted:
            self.face_weights = voro_info['face_areas']
            for i, f in enumerate(self.face_weights):
                self.face_weights[i] = [j/sum(f) for j in f]
        self.wigner_3j_cache = {}
        self.q_lm_cache = {}
        self.Q_lm_cache = {}

    def compute_q_lm(self, l, m):
        if (l, m) in self.q_lm_cache:
            return self.q_lm_cache[(l, m)]

        q_lm = np.zeros(self.n, dtype=np.complex128)
        for iatom in range(self.n):
            nn = len(self.neighbors[iatom])

            if nn == 0:
                q_lm[iatom] = np.nan + np.nan * 1j
            else:
                weights = 1.0/nn
                if self.weighted:
                    weights = np.asarray(self.face_weights[iatom])

                bvecs = self.neighbor_vectors[iatom]
                xy = bvecs[:, 0]**2 + bvecs[:, 1]**2
                theta = np.arctan2(bvecs[:, 1], bvecs[:, 0])
                phi = np.arctan2(np.sqrt(xy), bvecs[:, 2])
                t = sph_harm(m, l, theta, phi)

                q_lm[iatom] = np.sum(t*weights)

        self.q_lm_cache[(l, m)] = q_lm

        return q_lm


    def compute_Q_lm(self, l, m):
        if (l, m) in self.Q_lm_cache:
            return self.Q_lm_cache[(l, m)]

        q_lm = self.compute_q_lm(l, m)
        Q_lm = np.zeros_like(q_lm)
        for iatom in range(self.n):
            nn = len(self.neighbors[iatom])
            t = q_lm[iatom] + q_lm[self.neighbors[iatom]].sum()
            Q_lm[iatom] = t/(nn + 1)

        self.Q_lm_cache[(l, m)] = Q_lm

        return Q_lm

    def compute_q_l(self, l, coarse_grained=False):
        q_lm_func = self.compute_q_lm
        if coarse_grained:
            q_lm_func = self.compute_Q_lm

        q_l = 0
        for m in range(-l, l + 1):
            q_lm = q_lm_func(l, m)
            q_l += np.abs(q_lm) ** 2
        q_l = np.sqrt(4*np.pi/(2*l + 1)*q_l)

        return q_l

    def compute_w_l(self, l, coarse_grained=False):
        q_lm_func = self.compute_q_lm
        if coarse_grained:
            q_lm_func = self.compute_Q_lm

        if l not in self.wigner_3j_cache:
            w3j = _wigner_3j_values(l)
            self.wigner_3j_cache[l] = w3j
        else:
            w3j = self.wigner_3j_cache[l]

        w_l = 0
        for m1, m2, m3 in w3j.keys():
            q_lm1 = q_lm_func(l, m1)
            q_lm2 = q_lm_func(l, m2)
            q_lm3 = q_lm_func(l, m3)
            w = w3j[(m1, m2, m3)]
            w_l += w*q_lm1*q_lm2*q_lm3

        return w_l.real

    def compute_w_l_cap(self, l, coarse_grained=False):
        q_lm_func = self.compute_q_lm
        if coarse_grained:
            q_lm_func = self.compute_Q_lm

        q_l = 0
        for m in range(-l, l + 1):
            q_lm = q_lm_func(l, m)
            q_l += np.abs(q_lm)**2
        w_l = self.compute_w_l(l, coarse_grained=coarse_grained)
        w_l_cap = w_l*np.power(q_l, -3/2)

        return w_l_cap


if __name__ == "__main__":
    from ase.build import bulk
    from voroana import voronoi_analysis

    a = bulk("Fe")*(2, 2, 2)
    boo = BondOrientOrder(a, voronoi_analysis(a, outputs='snN'))
    print(boo.compute_q_l(4)[0])
    print(boo.compute_q_l(6)[0])
