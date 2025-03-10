import numpy as np
from world import World


class Body:

    def __init__(self, is_dynamic=True):
        super().__init__()

        self.density = 1.0

        # values expresses the affine mass
        self.M = None
        self.gravity = np.array([0, -9.81])

        self.volume = 0  # volume of the body

        self.kappa = 100000  # 10MPa

        self.q       = np.array([0., 0., 1., 0., 0., 1.])        # row vec
        self.qdot    = np.array([0., 0., 0., 0., 0., 0.])     # row vec
        self.q_tilde = np.array([0., 0., 1., 0., 0., 1.])  # row vec

        self.vertices = None
        self.current_vertices = None

        self.world = None

        self.is_dynamic = is_dynamic


    def set_vertices(self, vertices):
        self.vertices = vertices
        self.current_vertices = np.copy(vertices)
        self.update_mass_properties()

    def read_body_from_json(self, json_file):
        import json
        with open(json_file) as f:
            data = json.load(f)
            self.density = data['density']
            self.vertices = np.array(data['vertices'])
            self.kappa = data['kappa']

        self.current_vertices = np.copy(self.vertices)
        self.update_mass_properties()


    def update_mass_properties(self):
        from common import calculate_mass_properties
        self.volume, M_1, M_x, M_y, M_x2, M_y2, M_xy = calculate_mass_properties(self.vertices)
        self.M = np.array([[M_1, 0,   M_x,  M_y,  0 ,      0],
                           [0,   M_1, 0,    0,    M_y,   M_x],
                           [M_x, 0,   M_x2, M_xy, 0,      0],
                           [M_y, 0,   M_xy, M_y2, 0,       0],
                           [0,   M_x, 0,    0,    M_x2, M_xy],
                           [0,   M_y, 0,    0,    M_xy, M_y2]])

    def update_q_tilde(self):
        affine_gravity = np.array([self.gravity[0], self.gravity[1], 0, 0, 0, 0])
        self.q_tilde = self.q + self.world.delta_t * self.qdot + self.world.delta_t ** 2 * affine_gravity

    def update_vertices(self):

        R = np.array([[self.q[2], self.q[3]], [self.q[4], self.q[5]]])
        p = np.array([self.q[0], self.q[1]])
        debug = True
        for i in range(len(self.vertices)):
            self.current_vertices[i] = np.dot(R, self.vertices[i]) + p

    def set_position(self, x, y):
        self.q[0] = x
        self.q[1] = y


class Box(Body):

    def __init__(self, is_dynamic=True):
        super().__init__(is_dynamic)

    def set_as_box(self, width, height):

        hf_w = width / 2
        hf_h = height / 2
        self.set_vertices(np.array([[-hf_w, hf_h], [-hf_w, -hf_h], [hf_w, -hf_h], [hf_w, hf_h]]))

