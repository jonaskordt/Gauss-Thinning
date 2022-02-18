import math
import numpy as np
from numpy import linalg as LA
import igl
from gauss_thinning import *


def expand_path(path, r, V, F):
    nv = len(V)
    adj = triangle_adjacency(F, nv)
    B = igl.barycenter(V, F)
    stack = []
    flag = [-1 for _ in range(len(B))]
    flag_result = [0 for _ in range(len(B))]
    result = []

    for i in path:
        stack.append(i)
        flag[i] = i

        while len(stack) > 0:
            id = stack.pop()

            if flag_result[id] == 0:
                result.append(id)
                flag_result[id] = 1

            for j in adj[id]:
                if flag[j] != i and LA.norm(B[i] - B[j]) < r:
                    stack.append(j)
                    flag[j] = i

    return result


def collect_active_neighbors(adj, B, active_triangles, N, r, nr):
    stack = []
    flag = [-1 for _ in range(len(B))]
    result = [[] for _ in range(len(B))]
    normal_cone_threshold = math.cos(nr * math.pi / 180)

    for i in active_triangles:
        stack.append(i)
        flag[i] = i

        while len(stack) > 0:
            id = stack.pop()

            result[i].append(id)

            for j in adj[id]:
                if (
                    flag[j] != i
                    and LA.norm(B[i] - B[j]) < r
                    and np.dot(N[i], N[j]) > normal_cone_threshold
                ):
                    stack.append(j)
                    flag[j] = i

    return result


def fit_active_normals(nbh, N, cosine_threshold, sigma=1):
    nf = len(nbh)
    N2 = np.zeros((nf, 3))
    angle_threshold = cosine_threshold * math.pi / 180

    for i in range(nf):
        nbi = nbh[i]

        if len(nbi) == 0:
            N2[i] = N[i]
            continue

        NN = np.zeros((len(nbi), 3))

        for k in range(len(nbi)):
            NN[k] = N[nbi[k]]

        W = -np.eye(len(nbi))

        if sigma < 10:
            for j in range(len(nbi)):
                dot = np.dot(NN[0], NN[j])
                if dot >= 1:
                    W[j][j] = 1
                elif dot < 0:
                    W[j][j] = 0
                else:
                    W[j][j] = math.exp(
                        -math.pow(math.acos(dot) / angle_threshold / sigma, 2)
                    )

        else:
            W = np.eye(len(nbi))

        _, _, frame = LA.svd(np.dot(NN.T, np.dot(W, NN)))
        left_cols = frame[0:2, :].T
        N2[i] = np.dot(np.dot(left_cols, left_cols.T), N[i].T)
        norm = LA.norm(N2[i])
        if norm > 0:
            N2[i] = N2[i] / LA.norm(N2[i])
        else:
            N2[i] = N[i]

    return N2


async def constraint_gauss_thinning(
    V,
    F,
    active_triangles,
    num_iterations=100,
    min_cone_angle=2.5,
    smooth=1e-5,
    start_angle=25,
    radius=0.1,
    sigma=2,
    callback=None,
):
    cone_angle = start_angle
    eps = 1e-3

    nv = len(V)
    center(V)

    TT = triangle_adjacency(F, nv)
    C = igl.cotmatrix_entries(V, F)
    L = igl.cotmatrix(V, F)
    M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_BARYCENTRIC).todense()
    B = igl.barycenter(V, F)

    if smooth:
        A = -L + smooth * np.dot(L.T, L) + eps * M
    else:
        A = -L + eps * M

    V_dash = np.copy(V)

    for i in range(num_iterations):
        N = igl.per_face_normals(V_dash, F, np.array([1, 0, 0], dtype=float))
        B = igl.barycenter(V_dash, F)

        nbhs = collect_active_neighbors(TT, B, active_triangles, N, radius, cone_angle)
        if cone_angle > min_cone_angle:
            cone_angle *= 0.95

        N2 = fit_active_normals(nbhs, N, cone_angle, sigma)
        rot = find_rotations(N, N2)
        b = assemble_RHS(C, V_dash, F, rot)

        V_dash = np.asarray(LA.solve(A, eps * np.dot(M, V_dash) - b))

        if callback is not None:
            await callback(V_dash)

        print(f"Finished iteration {i}.")

    return V_dash


async def local_gauss_thinning(
    V,
    F,
    path,
    brush_size=0.05,
    **kwargs,
):
    active_triangles = expand_path(path, brush_size, V, F)
    return (
        await constraint_gauss_thinning(V, F, active_triangles, **kwargs),
        active_triangles,
    )
