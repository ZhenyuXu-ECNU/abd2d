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


class PointEdgeDistanceType(Enum):
    P_E0 = 1
    P_E1 = 2
    P_E = 3
    AUTO = 4


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


def point_edge_distance(p, e0, e1):

    distance_type = get_point_edge_distance_type(p, e0, e1)

    match distance_type:
        case PointEdgeDistanceType.P_E0:
            return point_point_distance(p, e0)
        case PointEdgeDistanceType.P_E1:
            return point_point_distance(p, e1)
        case PointEdgeDistanceType.P_E:
            return point_line_distance(p, e0, e1)
        case PointEdgeDistanceType.AUTO:
            raise ValueError("Invalid distance type of pint-edge!")

