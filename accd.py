# from cmath import sqrt
import sys
import numpy as np
from distance import point_edge_distance, edge_edge_distance
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


    def edge_edge_ccd(self, ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, min_distance, tmax):

        initial_distance = edge_edge_distance(ea0_t0, ea1_t0, eb0_t0, eb1_t0)

        if initial_distance <= min_distance * min_distance:
            warnings.warn(f"Initial distance {initial_distance} ≤ d_min={min_distance} is less than min_distance!")
            return [True, 0]

        dea0 = ea0_t1 - ea0_t0
        dea1 = ea1_t1 - ea1_t0
        deb0 = eb0_t1 - eb0_t0
        deb1 = eb1_t1 - eb1_t0
        subtract_mean(dea0, dea1, deb0, deb1)

        max_disp_mag = sqrt(max(dea0.dot(dea0), dea1.dot(dea1))) + sqrt(max(deb0.dot(deb0), deb1.dot(deb1)))

        if max_disp_mag == 0:
            return [False, None]

        min_distance_sq = min_distance * min_distance

        def distance_squared(x):
            d_sq = edge_edge_distance(x[0:2], x[2:4], x[4:6], x[6:8])
            if d_sq - min_distance_sq <= 0:
                d_sq = min(np.square(x[0:2] - x[4:6]), np.square(x[0:2] - x[6:8]),
                           np.square(x[2:4] - x[4:6]), np.square(x[2:4] - x[6:8]))
            return d_sq

        x = np.concatenate((ea0_t0, ea1_t0, eb0_t0, eb1_t0))
        dx = np.concatenate((dea0, dea1, deb0, deb1))

        return self.additive_ccd(x, dx, distance_squared, max_disp_mag, min_distance, tmax)