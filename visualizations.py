import numpy as np
import igl
from numpy import linalg as LA
import math


def angle_sum_colors(v, f, epsilon=1):
    vf, ni = igl.vertex_triangle_adjacency(f, v.shape[0])
    angle_sums = np.zeros(v.shape[0])
    for i in range(v.shape[0]):
        faces = f[vf[ni[i] : ni[i + 1]]]
        angle_sum = 0
        vertex = v[i]
        for face in faces:
            sides = [v[point] - vertex for point in face if point != i]
            norm_sides = [side / LA.norm(side) for side in sides]
            angle_sum = angle_sum + math.acos(np.dot(norm_sides[0], norm_sides[1]))
        angle_sums[i] = angle_sum / (math.pi) * 180

    return np.array(
        [
            0 if a < 360 - epsilon else 1 if a > 360 + epsilon else 0.5
            for a in angle_sums
        ]
    )


def developability_colors(v, f):
    nf = igl.per_face_normals(v, f, np.array([0.0, 0.0, 0.0]))
    vf, ni = igl.vertex_triangle_adjacency(f, v.shape[0])
    developability = np.zeros(v.shape[0])
    for i in range(0, v.shape[0]):
        dmat = nf[vf[ni[i] : ni[i + 1]]]
        developability[i] = min(np.linalg.eigvals(np.dot(dmat.T, dmat)))
    return developability
