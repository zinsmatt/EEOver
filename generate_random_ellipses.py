#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 23:11:00 2020

@author: mzins
"""

from ellcv.types import Ellipse
import numpy as np


out_file = "random_input.txt"
pos_range = [-7, 7]
axes_range = [1, 5]


def random_ellipse():
    ax = np.random.uniform(axes_range[0], axes_range[1], size=(2))
    ang = np.random.uniform(-np.pi, np.pi)
    pos = np.random.uniform(pos_range[0], pos_range[1], size=(2))
    return Ellipse.compose(ax, ang, pos)

def to_string(ell):
    ax, ang, pos = ell.decompose()
    return "%f %f %f %f %f" % (ax[0], ax[1], pos[0], pos[1], ang)

N = 100


lines = []
for ti in range(N):
    scale = np.random.uniform(0, 100, size=(1))
    ell0 = random_ellipse().full_scale(scale)
    ell1 = random_ellipse().full_scale(scale)

    line = " ".join([str(ti), to_string(ell0), to_string(ell1)]) + "\n"
    lines.append(line)

with open(out_file, "w") as fout:
    fout.writelines(lines)


