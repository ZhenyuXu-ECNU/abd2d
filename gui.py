import pygame
import numpy as np


class Gui:

    def __init__(self, resolution):
        self.screen = pygame.display.set_mode(resolution)
        self.resolution = np.array(resolution)
        self.offset = self.resolution / 2
        self.scale = 200
        self.running = True

    def gui_event_handler(self, event):
        if event.type == pygame.MOUSEWHEEL:
            if event.y > 0:
                self.scale *= 1.1
            else:
                self.scale *= 0.9
            return self.running

        if event.type == pygame.MOUSEMOTION:
            lef, mid, right = pygame.mouse.get_pressed()
            if right:
                self.offset[0] += event.rel[0]
                self.offset[1] -= event.rel[1]

    def screen_projection(self, x):
        return [self.offset[0] + self.scale * x[0], self.resolution[1] - (self.offset[1] + self.scale * x[1])]

    def draw_body(self, body):
        if body.vertices is None:
            return

        for i in range(len(body.vertices)):
            pygame.draw.aaline(self.screen, (0, 0, 255), self.screen_projection(body.vertices[i]),
                               self.screen_projection(body.vertices[(i + 1) % len(body.vertices)]))

