import numpy as np

def calc_inertia(m, ignore_hydrogen=False):
    arr = np.zeros([3, 3])
    for aa in m.Atoms:
        if ignore_hydrogen and aa.elem == "H":
            continue
        arr[0][0] += aa.mass * (aa.coor[1] ** 2 + aa.coor[2] ** 2)
        arr[0][1] -= aa.mass * aa.coor[0] * aa.coor[1]
        arr[0][2] -= aa.mass * aa.coor[0] * aa.coor[2]
        arr[1][0] -= aa.mass * aa.coor[1] * aa.coor[0]
        arr[1][1] += aa.mass * (aa.coor[0] ** 2 + aa.coor[2] ** 2)
        arr[1][2] -= aa.mass * aa.coor[1] * aa.coor[2]
        arr[2][0] -= aa.mass * aa.coor[2] * aa.coor[0]
        arr[2][1] -= aa.mass * aa.coor[2] * aa.coor[1]
        arr[2][2] += aa.mass * (aa.coor[0] ** 2 + aa.coor[1] ** 2)
    tt, qq = np.linalg.eig(arr)
    mod = np.linalg.norm(tt, ord=2)
    return (tt, mod)
