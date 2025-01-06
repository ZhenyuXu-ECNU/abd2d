import numpy as np


class Body:

    def __init__(self):
        self.density = 1.0

        # values expresses the affine mass
        self.M_1 = 0
        self.M_x = 0
        self.M_y = 0
        self.M_x2 = 0
        self.M_y2 = 0
        self.M_xy = 0

        self.volume = 0  # volume of the body

        self.kappa = 100000  # 10MPa

        self.q = np.array([0, 0, 1, 0, 0, 1])
        self.qdot = np.array([0, 0, 0, 0, 0, 0])

        self.vertices = None

    def read_body_from_json(self, json_file):
        import json
        with open(json_file) as f:
            data = json.load(f)
            self.density = data['density']
            self.vertices = np.array(data['vertices'])
            self.kappa = data['kappa']

        from common import calculate_mass_properties
        self.volume, self.M_1, self.M_x, self.M_y, self.M_x2, self.M_y2, self.M_xy = calculate_mass_properties(self.vertices)

    def orthogonality_potential(self):

        a11 = self.q[2]
        a12 = self.q[3]
        a21 = self.q[4]
        a22 = self.q[5]

        a1 = np.matrix([a11, a12]).transpose()
        a2 = np.matrix([a21, a22]).transpose()

        coeff = self.volume * self.kappa

        v = coeff * (((a1.transpose() * a1).item() - 1)**2 + ((a2.transpose() * a2).item() - 1)**2)
        v += 2 * coeff * (a1.transpose() * a2).item()**2

        dvda1 = 2 * coeff * (2 * ((a1.transpose() * a1).item() - 1) * a1 + 2 * (a2 * a2.transpose()) * a1)
        dvda2 = 2 * coeff * (2 * ((a2.transpose() * a2).item() - 1) * a2 + 2 * (a1 * a1.transpose()) * a2)

        hvda1da1 = 2 * coeff * (4 * a1 * a1.transpose() + 2 * ((a1.transpose() * a1).item() - 1) * np.identity(2) + 2 * (a2 * a2.transpose()))
        hvda2da2 = 2 * coeff * (4 * a2 * a2.transpose() + 2 * ((a2.transpose() * a2).item() - 1) * np.identity(2) + 2 * (a1 * a1.transpose()))

        hvda1da2 = 2 * coeff * np.matrix([[4 * a11 * a21 + 2 * a12 * a22, 2 * a12 * a21],
                                          [2 * a11 * a22, 2 * a11 * a21 + 4 * a12 * a22]])

        return v, dvda1, dvda2, hvda1da1, hvda2da2, hvda1da2




