import math
import numpy as np

def normal_vector_plane(point1, point2, point3):
    # 计算向量 v1 和 v2
    point1 = np.array(point1)
    point2 = np.array(point2)
    point3 = np.array(point3)
    v1 = point2 - point1
    v2 = point3 - point1
    #v1 = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
    #v2 = [point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2]]
    
    # 计算叉积
    #cross_product = [
    #    v1[1]*v2[2] - v1[2]*v2[1],
    #    v1[2]*v2[0] - v1[0]*v2[2],
    #    v1[0]*v2[1] - v1[1]*v2[0]
    #]
    cross_product = np.cross(v1,v2)
    length = np.linalg.norm(cross_product)
    normal_vector = cross_product / length

    # 计算长度
    #length = math.sqrt(cross_product[0]**2 + cross_product[1]**2 + cross_product[2]**2)
    
    # 归一化向量
    #normal_vector = [cross_product[0]/length, cross_product[1]/length, cross_product[2]/length]
    
    return normal_vector

def plane_to_plane_angle(plane1,plane2,normal_vector=False,cross_line = False):
    n1 = normal_vector_plane(*plane1)
    n2 = normal_vector_plane(*plane2)
    dot_product = np.dot(n1,n2)
    length1 = np.linalg.norm(n1)
    length2 = np.linalg.norm(n2)
    cosine = dot_product / length1 / length2
    if cosine > 1.0:
        cosine = 1.0
    elif cosine < -1.0:
        cosine = -1.0
    angle = math.degrees(math.acos(cosine))
    if cross_line:
        cl = np.cross(n1,n2)
        if normal_vector:
            return angle, cl, n1, n2
        else:
            return angle, cl
    else:
        if normal_vector:
            return angle, n1, n2
        else:
            return angle

def calculate_dihedral(a,b,c,d):
    angle, n1, n2 = plane_to_plane_angle([a,b,c],[b,c,d],normal_vector=True)
    rad = np.array(d) - np.array(a)
    if np.dot(n1,rad) < 0.0:
        return angle * -1
    return angle

def rotation_dihedral(b, c, d, tht, radian=False):
    """
    Rotate the atom d along the axis b-c by magnitude of angle. The rotation is
    anti-clockwise if angle > 0.

    Args:
        b: [float, float, float], cartesian coordinate of the atom in rotation
            axis
        c: [float, float, float], cartesian coordinate of the atom in rotation
            axis
        d: [float, float, float], cartesian coordinate of the atom to be rotated
        tht (float): value of angle to be changed, default unit: degrees
        radian (boolean): True when the unit of input angle is in radial

    Returns:
        [float, float, float], rotated coordiante of atom d
    """
    bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]
    r = (bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]) ** 0.5
    nbc = [bc[0] / r, bc[1] / r, bc[2] / r]
    if not radian:
        tht = tht / 180 * np.pi
    a1 = math.cos(tht)
    a2 = 1 - a1
    a3 = math.sin(tht)
    nx = nbc[0]
    ny = nbc[1]
    nz = nbc[2]
    dx = d[0] - b[0]
    dy = d[1] - b[1]
    dz = d[2] - b[2]
    m = [
        [a1 + a2 * nx * nx, a2 * nx * ny - a3 * nz, a2 * nx * nz + a3 * ny],
        [a2 * ny * nx + a3 * nz, a1 + a2 * ny * ny, a2 * ny * nz - a3 * nx],
        [a2 * nz * nx - a3 * ny, a2 * nz * ny + a3 * nx, a1 + a2 * nz * nz],
    ]
    ddx = m[0][0] * dx + m[0][1] * dy + m[0][2] * dz + b[0]
    ddy = m[1][0] * dx + m[1][1] * dy + m[1][2] * dz + b[1]
    ddz = m[2][0] * dx + m[2][1] * dy + m[2][2] * dz + b[2]
    dd = [ddx, ddy, ddz]
    return dd
