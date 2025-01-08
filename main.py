import copy

import pygame
import numpy as np
from body import Body
from world import World
from ipc_energy import TotalIpc, OrthogonalEnergy, InertiaEnergy
import numpy.linalg as LA

from gui import Gui
from cmath import inf


tol = 1e-3


if __name__ == '__main__':

    world = World()

    box = Body()
    box.read_body_from_json('fixtures/box.json')

    world.add_body(box)

    o1 = OrthogonalEnergy.val(box)
    o2 = OrthogonalEnergy.grad(box)  # row vec
    o3 = OrthogonalEnergy.hess(box)  # row major

    i1 = InertiaEnergy.val(box)
    i2 = InertiaEnergy.grad(box)
    i3 = InertiaEnergy.hess(box)

    gui = Gui([1920, 1080])
    running = True

    while running:

        for event in pygame.event.get():
            # check for closing window
            if event.type == pygame.QUIT:
                running = False

            gui.gui_event_handler(event)  # check for other events

        gui.screen.fill((255, 255, 255))

        box.update_vertices()

        gui.draw_body(box)
        pygame.display.flip()

        # ipc solve
        box.update_q_tilde()

        E_last = TotalIpc.val(world.bodies)
        p_box = TotalIpc.search_dir(box)

        q_last = copy.deepcopy(box.q)

        while LA.norm(p_box, inf) / world.delta_t > tol:

            # line search
            alpha = 1
            while TotalIpc.val(world.bodies) > E_last:
                alpha /= 2

            box.q += alpha * p_box

            E_last = TotalIpc.val(world.bodies)
            p_box = TotalIpc.search_dir(box)

        print('box q:', box.q)
        box.qdot = (box.q - q_last) / world.delta_t

        pygame.time.wait(int(world.delta_t * 1000))
