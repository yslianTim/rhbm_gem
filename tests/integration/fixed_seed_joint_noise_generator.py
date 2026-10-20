"""Deterministic map and noise generators retained from the closed study."""
import math

import numpy as np


SEEDS = tuple(range(20260921, 20260941))
AXIS = np.arange(25) * .5 - 6
Z, Y, X = np.meshgrid(AXIS, AXIS, AXIS, indexing='ij')
TARGET = (X * X + Y * Y + Z * Z <= 6.25).ravel()


def clean_map(distance, shift):
    values = np.zeros(X.shape)
    for x, a, b, c in [(0, 2, .5, .2), (distance + shift, 2.3, .55, .15)]:
        square = (X - x) ** 2 + Y ** 2 + Z ** 2
        radius = np.sqrt(square)
        keep = square <= 6.25
        gaussian = np.exp(-square / (2 * b * b)) / (2 * math.pi * b * b) ** 1.5
        charge = np.zeros(X.shape)
        center = radius < 1e-5
        charge[center] = math.sqrt(2 / math.pi) / b
        nonzero = ~center
        charge[nonzero] = np.fromiter(
            (math.erf(float(t) / (math.sqrt(2) * b)) / float(t) for t in radius[nonzero]), float)
        values[keep] += (a * gaussian + c * charge)[keep]
    return values.ravel()


def noise_field(seed, kind):
    rng = np.random.Generator(np.random.PCG64(seed))
    if kind == 'iid':
        return rng.normal(size=(25, 25, 25)).ravel()
    field = rng.normal(size=(31, 31, 31))
    kernel = np.exp(-.5 * np.arange(-3, 4, dtype=float) ** 2)
    kernel /= np.sqrt(np.sum(kernel * kernel))
    for axis in range(3):
        field = np.apply_along_axis(lambda row: np.convolve(row, kernel, mode='valid'), axis, field)
    return field.ravel()
