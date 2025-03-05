import numpy as np
from enum import Enum
import warnings


# using squared form distance
def point_point_distance(p0, p1):
    return np.dot(p0 - p1, p0 - p1)


# using squared form distance
def point_line_distance(p, e0, e1):
    e = e1 - e0
    numerator = e[1] * p[0] - e[0] * p[1] + e1[0] * e0[1] - e1[1] * e0[0]
    return numerator * numerator / np.dot(e, e)


def line_line_distance(ea0, ea1, eb0, eb1):
    n = (ea1 - ea0).cross(eb1 - eb0)
    line_to_line = (eb0 - ea0).dot(n)
    return line_to_line * line_to_line / n.dot(n)


class PointEdgeDistanceType(Enum):
    P_E0 = 1
    P_E1 = 2
    P_E = 3
    AUTO = 4


class EdgeEdgeDistanceType(Enum):
    EA0_EB0 = 1,
    EA0_EB1 = 2,
    EA1_EB0 = 3,
    EA1_EB1 = 4,
    EA_EB0 = 5,
    EA_EB1 = 6,
    EA0_EB = 7,
    EA1_EB = 8,
    EA_EB = 9,
    AUTO = 10


def get_point_edge_distance_type(p, e0, e1):

    e = e1 - e0  # the edge vector

    e_length_sqr = np.dot(e, e)

    if e_length_sqr == 0:
        warnings.warn("Degenerate edge in point_edge_distance_type!", UserWarning)
        return PointEdgeDistanceType.P_E0  # arbitrary choice

    ratio = np.dot(p - e0, e) / e_length_sqr

    if ratio < 0:
        return PointEdgeDistanceType.P_E0
    elif ratio > 1:
        return PointEdgeDistanceType.P_E1
    else:
        return PointEdgeDistanceType.P_E


def get_edge_edge_distance_type(ea0, ea1, eb0, eb1):

    parallel_threshold = 1e-20

    u = ea1 - ea0
    v = eb1 - eb0
    w = ea0 - eb0

    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w)
    e = np.dot(v, w)
    D = a * c - b * b

    # degenerated case
    if a == 0 and c == 0:
        return EdgeEdgeDistanceType.EA0_EB0
    elif a == 0:
        return EdgeEdgeDistanceType.EA0_EB
    elif c == 0:
        return EdgeEdgeDistanceType.EA_EB0

    parallel_tolerance = parallel_threshold * max(1.0, a * c)

    if np.square(u.cross(v)) < parallel_tolerance:
        return edge_edge_parallel_distance_type(ea0, ea1, eb0, eb1)

    default_case = EdgeEdgeDistanceType.EA_EB

    s_N = b * e - c * d
    if s_N <= 0:
        t_N = e
        t_D = c
        default_case = EdgeEdgeDistanceType.EA0_EB0
    elif s_N >= D:
        t_N = e + b
        t_D = c
        default_case = EdgeEdgeDistanceType.EA1_EB
    else:
        t_N = a * e - b * d
        t_D = D
        if 0 < t_N < t_D and np.square(u.cross(v) < parallel_tolerance):
            if s_N < D/2:
                t_N = e
                t_D = c
                default_case = EdgeEdgeDistanceType.EA0_EB
            else:
                t_N = e + b
                t_D = c
                default_case = EdgeEdgeDistanceType.EA1_EB

    if t_N <= 0:
        if -d <= 0:
            return EdgeEdgeDistanceType.EA0_EB0
        elif -d >= a:
            return EdgeEdgeDistanceType.EA1_EB0
        else:
            return EdgeEdgeDistanceType.EA_EB0
    elif t_N >= t_D:
        if -d + b <= 0:
            return EdgeEdgeDistanceType.EA0_EB1
        elif -d + b >= a:
            return EdgeEdgeDistanceType.EA1_EB1
        else:
            return EdgeEdgeDistanceType.EA_EB1
    return default_case


def point_edge_distance(p, e0, e1):

    distance_type = get_point_edge_distance_type(p, e0, e1)

    match distance_type:
        case PointEdgeDistanceType.P_E0:
            return point_point_distance(p, e0)
        case PointEdgeDistanceType.P_E1:
            return point_point_distance(p, e1)
        case PointEdgeDistanceType.P_E:
            return point_line_distance(p, e0, e1)
        case _:
            raise ValueError("Invalid distance type of pint-edge!")


def edge_edge_parallel_distance_type(ea0, ea1, eb0, eb1):
    ea = ea1 - ea0
    alpha = (eb0 - ea0).dot(ea) / ea.dot(ea)
    beta = (eb1 - ea0).dot(ea) / ea.dot(ea)

    if alpha < 0:
        eac =  2 if 0 <= beta &beta <= 1 else 0
        ebc = 0  if beta <= alpha else 1 if beta <= 1 else 2
    elif alpha > 1:
        eac = 2 if 0 <= beta & beta <= 1 else 1
        ebc = 0 if beta >= alpha else 1 if 0 <= beta else 2
    else:
        eac = 2
        ebc = 0
    return EdgeEdgeDistanceType((eac << 1 or ebc) if ebc < 2 else (6 + eac))


def edge_edge_distance(ea0, ea1, eb0, eb1):

    distance_type = get_edge_edge_distance_type(ea0, ea1, eb0, eb1)
    match distance_type:
        case EdgeEdgeDistanceType.EA0_EB0:
            return point_point_distance(ea0, eb0)
        case EdgeEdgeDistanceType.EA0_EB1:
            return point_point_distance(ea0, eb1)
        case EdgeEdgeDistanceType.EA1_EB0:
            return point_point_distance(ea1, eb0)
        case EdgeEdgeDistanceType.EA1_EB1:
            return point_point_distance(ea1, eb1)
        case EdgeEdgeDistanceType.EA_EB0:
            return point_line_distance(eb0, ea0, ea1)
        case EdgeEdgeDistanceType.EA_EB1:
            return point_line_distance(eb1, ea0, ea1)
        case EdgeEdgeDistanceType.EA0_EB:
            return point_line_distance(ea0, eb0, eb1)
        case EdgeEdgeDistanceType.EA1_EB:
            return point_line_distance(ea1, eb0, eb1)
        case EdgeEdgeDistanceType.EA_EB:
            return line_line_distance(ea0, ea1, eb0, eb1)
        case _:
            raise ValueError("Invalid distance type of edge-edge!")

