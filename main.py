import copy

import pygame
import numpy as np
from body import Body, Box
from world import World
from ipc_energy import TotalIpc, OrthogonalEnergy, InertiaEnergy
import numpy.linalg as LA

from gui import Gui
from cmath import inf


tol = 1e-3


if __name__ == '__main__':

    world = World()

    box = Box()
    box.set_as_box(1.0, 1.0)
    # box.read_body_from_json('fixtures/box.json')

    background = Box()
    background.set_as_box(60.0, 0.5)
    background.is_dynamic = False
    background.set_position(0, -2.0)

    world.add_body(box)
    world.add_body(background)

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
        background.update_vertices()

        gui.draw_body(box)
        gui.draw_body(background)
        pygame.display.flip()

        ######################################## ipc solve ########################################
        box.update_q_tilde()
        background.update_q_tilde()

        E_last = TotalIpc.val(world.bodies)
        p_box = TotalIpc.search_dir(box)
        p_background = TotalIpc.search_dir(background)

        q_last_box = copy.deepcopy(box.q)
        q_last_background = copy.deepcopy(background.q)

        while LA.norm(np.concatenate((p_box, p_background)), inf) / world.delta_t > tol:

            # line search
            alpha = 1
            while TotalIpc.val(world.bodies) > E_last:
                alpha /= 2

            box.q += alpha * p_box
            background.q += alpha * p_background

            E_last = TotalIpc.val(world.bodies)
            p_box = TotalIpc.search_dir(box)
            p_background = TotalIpc.search_dir(background)

            debug = 1

        box.qdot = (box.q - q_last_box) / world.delta_t
        print('box.qdot =', box.qdot)
        background.qdot = (background.q - q_last_background) / world.delta_t

        pygame.time.wait(int(world.delta_t * 1000))
