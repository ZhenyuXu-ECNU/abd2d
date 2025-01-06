

def calculate_mass_properties(vertices):
    volumn = 0
    M_1 = 0
    M_x = 0
    M_y = 0
    M_x2 = 0
    M_y2 = 0
    M_xy = 0

    for i in range(len(vertices)):
        x0 = vertices[i][0]
        y0 = vertices[i][1]

        x1 = vertices[(i + 1) % len(vertices)][0]
        y1 = vertices[(i + 1) % len(vertices)][1]

        nx = y1 - y0
        ny = x0 - x1

        t0 = x0 + x1
        t1 = x0 * x0 + x0 * x1 + x1 * x1
        t2 = y0 * y0 + y0 * y1 + y1 * y1
        t3 = (x0 * x0 + x1 * x1) * (x0 + x1)
        t4 = (y0 * y0 + y1 * y1) * (y0 + y1)
        t5 = x0 * x0 * (3 * y0 + y1) + x0 * x1 * (2 * y0 + 2 * y1) + x1 * x1 * (y0 + 3 * y1)

        M_1 += (1 / 2) * nx * t0
        M_x += (1 / 6) * nx * t1
        M_y += (1 / 6) * ny * t2
        M_x2 += (1 / 12) * nx * t3
        M_y2 += (1 / 12) * ny * t4
        M_xy += (1 / 24) * nx * t5

    return M_1, M_1, M_x, M_y, M_x2, M_y2, M_xy