import numpy as np
import pymorton


class BVH:

    def __init__(self):
        self.boxlist = []
        self.new2old = []
        self.n_corners = -1

    def init(self, vertices, tol):

        corner_list = []
        for i in range(len(vertices)):
            v0 = vertices[i]
            v1 = vertices[(i + 1) % len(vertices)]

            min_corner = np.minimum(v0, v1) - tol
            max_corner = np.maximum(v0, v1) + tol

            corner_list.append((min_corner, max_corner))

        self.init(corner_list)

    def init(self, corner_list):
        n_corners = len(corner_list)

        box_centers = np.array([0.5 * (corner_list[i][0] + corner_list[i][1]) for i in range(n_corners)])
        min_boxes = np.min(box_centers, axis=0)
        max_boxes = np.max(box_centers, axis=0)

        center = 0.5 * (min_boxes + max_boxes)

        box_centers = box_centers - center

        scale_point = max_boxes - center

        # resize the boxes
        scale = np.linalg.norm(scale_point, ord=np.inf)
        if scale > 100:
            box_centers = box_centers / scale

        list = []
        multi = 1000

        for i in range(n_corners):
            list.append((i, pymorton.interleave2(1000 * box_centers[i].astype(int))))

        list = sorted(list, key=lambda x: x[1])

