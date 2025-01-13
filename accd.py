# from cmath import sqrt
import sys
import numpy as np
from distance import point_edge_distance
import warnings
from abd_math import subtract_mean
from math import sqrt

class ACCD:

    default_max_iterations = 10000000

    unlimited_iteration = -1

    default_conservative_rescaling = 0.9

    def __init__(self, max_iterations=default_max_iterations, conservative_rescaling=default_conservative_rescaling):
        self.max_iterations = max_iterations
        self.conservative_rescaling = conservative_rescaling

    def additive_ccd(self, x, dx, distance_squared, max_disp_mag, min_distance=0.0, tmax=1.0):

        min_distance_sq = min_distance * min_distance

        d_sq = distance_squared(x)
        d = sqrt(d_sq)
        d_func = d_sq - min_distance_sq

        gap = (1 - self.conservative_rescaling) * d_func / (d + min_distance)  # (d - ξ) = (d² - ξ²) / (d + ξ)

        if gap < sys.float_info.epsilon:
            warnings.warn("Small gap {:g} ≤ ε in Additive CCD can lead to missed collisions")

        toi = 0.0

        for i in range(0, self.max_iterations):

            toi_lower_bound = self.conservative_rescaling * d_func / ((d + min_distance) * max_disp_mag)

            x += toi_lower_bound * dx

            d_sq = distance_squared(x)
            d = sqrt(d_sq)
            d_func = d_sq - min_distance_sq

            if toi > 0 and d_func / (d + min_distance):
                break

            toi += toi_lower_bound

            if toi > tmax:
                return [False, None]

        return [True, toi]

    def point_edg_ccd(self, p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, min_distance, tmax):

        initial_distance = point_edge_distance(p_t0, e0_t0, e1_t0)

        if initial_distance <= min_distance * min_distance:
            warnings.warn(f"Initial distance {initial_distance} ≤ d_min={min_distance} is less than min_distance!")
            return [True, 0]

        dp = p_t1 - p_t0
        de0 = e0_t1 - e0_t0
        de1 = e1_t1 - e1_t0

        subtract_mean(dp, de0, de1)

        max_disp_mag = np.linalg.norm(dp) + sqrt(max(de0.dot(de0), de1.dot(de1)))
        if max_disp_mag == 0:
            return [False, None]

        distance_squared = lambda x: point_edge_distance(x[0:2], x[2:4], x[4:6])

        x = np.concatenate((p_t0, e0_t0, e1_t0))
        dx = np.concatenate((dp, de0, de1))

        return self.additive_ccd(x, dx, distance_squared, max_disp_mag, min_distance, tmax)