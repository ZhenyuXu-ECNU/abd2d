import numpy as np


def subtract_mean(*args):
    mean = np.mean(np.array(args), axis=0)

    for arg in args:
        arg -= mean

