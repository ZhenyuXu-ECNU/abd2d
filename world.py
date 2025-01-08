import numpy as np


class World:

    def __init__(self):
        self.gravity = np.array([0, -9.81])
        self.bodies = []
        self.delta_t = 0.01

    def add_body(self, body):
        self.bodies.append(body)
        body.gravity = self.gravity
        body.world = self
