import sys
import numpy as np
import numexpr as ne
from ase import Atoms
from voroana import voronoi_analysis
from voroana.io import read_dump, window_iter


class VoroFS:
    def __init__(self, dumpfile, voro_args={"outputs": "l"}):
        self.dumpfile = dumpfile
        self.voro_args = voro_args
        self.voro_cache = {}

    def voronoi(self, snapshot):
        if snapshot['index'] not in self.voro_cache:
            a = Atoms(positions=snapshot['positions'], cell=snapshot['box'], pbc=True)
            a.wrap()
            res = voronoi_analysis(a, **self.voro_args)
            self.voro_cache[snapshot['index']] = res['labels']
            return self.voro_cache[snapshot['index']]
        else:
            return self.voro_cache[snapshot['index']]

    def masker(self, snapshot, tag):
        if isinstance(tag, str) and tag == 'ALL':
            mask = np.array([True]*len(snapshot['types']))
            nmasked = mask.sum()
        else:
            if isinstance(tag, str):
                tagger = lambda x: x == tag
            labels = self.voronoi(snapshot)
            mask = np.array([tagger(i) for i in labels])
            nmasked = mask.sum()
        assert nmasked > 0

        return mask, nmasked

    def calc_sisf(self, q, dt, ncorr, nshift, tag="ALL", maxframes=0):
        dump = read_dump(self.dumpfile, maxframes=maxframes, dt=dt)
        ic = 0
        s = []
        for window in window_iter(dump, width=ncorr, stride=nshift):
            istart = window[0]['index']
            iend = window[-1]['index']
            sys.stdout.flush()
            mask, nmasked = self.masker(window[0], tag)
            print('VOROFS: {:3d}-th average, [{:5d} to {:5d}], n = {:5d}.'.format(ic, istart, iend, nmasked))
            res = np.empty((ncorr - 1, 2))
            res.fill(np.nan)
            for i in range(1, len(window)):
                dt = window[i]['time'] - window[0]['time']
                dpos = window[i]['positions'] - window[0]['positions']
                dpos = dpos[mask]
                sisf = ne.evaluate("cos(q*dpos)").mean(axis=1).sum()/nmasked
                res[i - 1, 0] = dt
                res[i - 1, 1] = sisf
            s.append(res)
            ic += 1
        s = np.stack(s, axis=2)
        res = np.zeros((ncorr - 1, 2))
        res[:, 0] = s[:, 0, 0]
        res[:, 1] = np.nanmean(s[:, 1, :], axis=1)

        return res

