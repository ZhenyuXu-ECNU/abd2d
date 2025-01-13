import pygame
import numpy as np


class Gui:

    def __init__(self, resolution):
        self.screen = pygame.display.set_mode(resolution)
        self.resolution = np.array(resolution)
        self.offset = self.resolution / 2
        self.scale = 100
        self.running = True
        self.point_size = 5

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

    def draw_body(self, body, color=(0, 0, 255)):
        if body.vertices is None:
            return

        for i in range(len(body.vertices)):
            pygame.draw.aaline(self.screen, color, self.screen_projection(body.current_vertices[i]),
                               self.screen_projection(body.current_vertices[(i + 1) % len(body.current_vertices)]))

    def draw_line(self, e1, e2, color=(0, 0, 0)):
        pygame.draw.aaline(self.screen, color, self.screen_projection(e1), self.screen_projection(e2))

    def draw_point(self, p, color=(0, 255, 255)):
        pygame.draw.circle(self.screen, color, self.screen_projection(p), self.point_size)
