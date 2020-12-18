from collections import deque
from itertools import islice
import sys
import numpy as np


def read_dump(fname, maxframes=0, dt=0.002):
    _ORTHO_BOX_FLAG = "ITEM: BOX BOUNDS pp pp pp"
    _DUMP_FLAGS = "ITEM: ATOMS id type xu yu zu"
    if maxframes == 0:
        maxframes = sys.maxsize
    iframe = 0
    f = open(fname, "r")
    tmp = [f.readline() for i in range(9)]
    # Check prerequisites
    if tmp[4].strip() != _ORTHO_BOX_FLAG:
        raise RuntimeError("Require periodic orthogonal box!")
    if tmp[8].strip() != _DUMP_FLAGS:
        raise RuntimeError("Require dump custom id type xu yu zu!")
    timestep_0 = int(tmp[1].strip())
    natoms_0 = int(tmp[3].strip())
    box_0 = np.array([list(map(float, i.strip().split())) for i in tmp[5:8]])
    box_0 = box_0[:, 1] -  box_0[:, 0]
    f.seek(0)
    def _read_one():
        if ('' == f.readline()):
            return None
        else:
            tmp = [f.readline() for i in range(8)]
            timestep = int(tmp[0].strip())
            natoms = int(tmp[2].strip())
            if not natoms == natoms_0:
                raise RuntimeError("Number of atoms changes during simulation!")
            box = np.array([list(map(float, i.strip().split())) for i in tmp[4:7]])
            box = box[:, 1] - box[:, 0]
            if not np.allclose(box, box_0):
                raise RuntimeError("Box changes during simulation!")
            types = np.zeros(natoms, dtype=np.int64)
            pos = np.zeros((natoms, 3), dtype=np.float64)
            for _ in range(natoms):
                l = f.readline()
                if '' == l:
                    raise RuntimeError("Unexpected file interruption!")
                l = l.strip().split()
                iid = int(l[0]) - 1
                types[iid] = int(l[1])
                pos[iid, 0] = float(l[2])
                pos[iid, 1] = float(l[3])
                pos[iid, 2] = float(l[4])

            res = {
                    "index"     : iframe,
                    "time"      : (timestep - timestep_0)*dt,
                    "box"       : box,
                    "types"     : types,
                    "positions" : pos
                  }
            return res

    _ = _read_one()
    iframe += 1
    while _ is not None and iframe <= maxframes:
        yield _
        _ = _read_one()
        iframe += 1


def consume(iterator, n=None):
    "Advance the iterator n-steps ahead. If n is None, consume entirely."
    # Use functions that consume iterators at C speed.
    # https://docs.python.org/zh-cn/3/library/itertools.html
    if n is None:
        # feed the entire iterator into a zero-length deque
        deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)


def window_iter(iterator, width=2, stride=1):
    assert width >= 1 and stride >= 1
    _it = iter(iterator)
    _window = None
    def one_window():
        nonlocal _it
        nonlocal _window
        if _window is None:
            _window = deque(islice(_it, width), width)
        else:
            if stride >= width:
                _window.clear()
                consume(_it, stride - width)
            else:
                for _ in range(min(stride, len(_window))):
                    _window.popleft()
            for f in islice(_it, min(stride, width)):
                _window.append(f)
        return _window

    _ = one_window()
    while len(_) > 1:
        yield _
        _ = one_window()


if __name__ == "__main__":
    for i in window_iter(range(13), 8, 3):
        print(i)
