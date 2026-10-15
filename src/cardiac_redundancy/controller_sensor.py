"""Causal time-weighted observation window; no future samples or padding."""
from collections import deque
import math
import numpy as np


class TrailingObservation:
    def __init__(self, window_s=.75):
        if not math.isfinite(window_s) or window_s <= 0:
            raise ValueError('Positive finite window required')
        self.window_s = window_s
        self.samples = deque()

    def update(self, t, pressure, native_flow):
        if any(not math.isfinite(v) for v in (t, pressure, native_flow)):
            raise ValueError('Finite observations required')
        if self.samples and t <= self.samples[-1][0]:
            raise ValueError('Strictly increasing sensor timestamps required')
        self.samples.append((t, pressure, native_flow))
        left = t-self.window_s
        # Keep the sample immediately before the window to interpolate its
        # boundary. Do not synthesize an unobserved prehistory at startup.
        while len(self.samples) > 2 and self.samples[1][0] <= left:
            self.samples.popleft()
        if self.samples[0][0] > left+1e-10:
            return None
        a = np.asarray(self.samples)
        times = np.r_[left, a[a[:, 0] > left, 0]]
        values = [np.interp(times, a[:, 0], a[:, j]) for j in (1, 2)]
        return tuple(float(np.trapezoid(v, times)/self.window_s) for v in values)
