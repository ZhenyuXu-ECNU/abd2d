import pygame
import numpy as np
from body import Body
from gui import Gui

if __name__ == '__main__':
    box = Body()
    box.read_body_from_json('fixtures/box.json')

    gui = Gui([1920, 1080])

    box.orthogonality_potential()
    running = True
    while running:

        for event in pygame.event.get():
            # check for closing window
            if event.type == pygame.QUIT:
                running = False

            gui.gui_event_handler(event)  # check for other events

        gui.screen.fill((255, 255, 255))

        gui.draw_body(box)

        pygame.display.flip()
