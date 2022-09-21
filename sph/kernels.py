from .const import *
import numpy as np

PI = np.pi


def is_close_enough(r_squared):
    return r_squared < DISTANCE_LIMIT_SQUARED


def poly_6(r):
    magic_number = 315.0 / 64.0 / PI / pow(DISTANCE_LIMIT, 9)
    return pow(DISTANCE_LIMIT - r * r, 3) * magic_number


def spiky(r):
    magig_number = 15.0 / PI / pow(DISTANCE_LIMIT, 6)
    return pow(DISTANCE_LIMIT - r, 3) * magig_number


def grad_spiky(r1, r2):
    dr = r1 - r2
    distance = np.linalg.norm(dr)
    if distance == 0:
        return np.array([0.0, 0.0, 0.0])
    magic_number = 45.0 / PI / pow(DISTANCE_LIMIT, 6)
    nominator = pow(DISTANCE_LIMIT - distance, 2) * magic_number
    return nominator * dr / distance


def viscosity(r):
    magic_number = 15.0 / 2.0 / PI / pow(DISTANCE_LIMIT, 3)
    ratio = r / DISTANCE_LIMIT
    nominator = -pow(ratio, 3) / 2 + pow(ratio, 2) + ratio / 2 - 1
    return magic_number * nominator


def laplacian_viscosity(r):
    magic_number = 45.0 / PI / pow(DISTANCE_LIMIT, 6)
    return magic_number * (DISTANCE_LIMIT - r)
