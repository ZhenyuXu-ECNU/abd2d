import pygame
from accd import ACCD
from gui import Gui
import numpy as np


delta_t = 0.1


if __name__ == '__main__':

    gui = Gui([1920, 1080])
    running = True

    point = np.array([-3.0, 0.0])
    e1, e2 = np.array([0, -1.0]), np.array([0, 1.0])

    accd = ACCD()

    is_collide_total = False

    while running:

        for event in pygame.event.get():
            # check for closing window
            if event.type == pygame.QUIT:
                running = False

            gui.gui_event_handler(event)  # check for other events

        gui.screen.fill((255, 255, 255))
        gui.draw_point(point)
        gui.draw_line(e1, e2)
        pygame.display.flip()

        ######################################## ccd solve ########################################

        ###########################################################################################
        pygame.time.wait(int(delta_t * 1000))

        is_collide, toi = accd.point_edg_ccd(point, e1, e2, point + delta_t * np.array([1.0, 0.0]), e1, e2, 0.0, 1.0)

        if is_collide:
            point += toi * delta_t * np.array([1.0, 0.0])
            delta_t = 0.0

        point += delta_t * np.array([1.0, 0.0])
