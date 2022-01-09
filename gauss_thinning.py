import igl
import meshplot as mp
import numpy as np
from numpy import linalg as LA
import math


def center(v):
    v -= np.mean(v, axis=0)
    v /= 2 * np.max(LA.norm(v, axis=1))


def triangle_adjacency(f, nv):
    vnbhs = [
        [] for _ in range(nv)
    ]  # for every vertex a list of faces the vertex is part of
    nf = len(f)

    for i in range(nf):
        for j in range(3):
            vnbhs[f[i][j]].append(i)

    flags = [-1 for _ in range(nf)]
    ret = [[] for _ in range(nf)]

    for i in range(nf):
        for j in range(3):
            for k in vnbhs[f[i, j]]:
                if k != i and flags[k] != i:
                    ret[i].append(k)
                    flags[k] = i

    return ret  # for every face a list of faces which are adjacent


def collectNeighbors(adj, V, N, r, nr):
    stack = []
    flag = [-1 for _ in range(len(V))]
    result = [[] for _ in range(len(V))]
    normal_cone_threshold = math.cos(nr * math.pi / 180)

    for i in range(len(V)):
        stack.append(i)
        flag[i] = i

        while len(stack) > 0:
            id = stack.pop()

            result[i].append(id)

            for j in adj[id]:
                if (
                    flag[j] != i
                    and LA.norm(V[i] - V[j]) < r
                    and np.dot(N[i], N[j]) > normal_cone_threshold
                ):
                    stack.append(j)
                    flag[j] = i

    return result


def fit_normals(nbh, N, cosine_threshold, sigma=1):
    nf = len(nbh)
    N2 = np.zeros((nf, 3))
    angle_threshold = cosine_threshold * math.pi / 180

    for i in range(nf):
        nbi = nbh[i]

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


def find_rotations(N0, N1):
    n = len(N0)
    rot = []

    for i in range(n):
        n1 = N0[i]
        n2 = N1[i]
        v = np.cross(n1, n2)
        c = np.dot(n1, n2)

        if c > -1 + 1e-8:
            coeff = 1 / (1 + c)
            v_x = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
            rot.append(np.eye(3) + v_x + coeff * np.dot(v_x, v_x))
        else:
            rot.append(-np.eye(3))

    return rot


def assemble_RHS(C, V, F, R):
    rhs = np.zeros((len(V), 3))

    for i in range(len(F)):
        for j in range(3):
            v0 = F[i][(j + 1) % 3]
            v1 = F[i][(j + 2) % 3]

            b = C[i][j] * np.dot(R[i], (V[v0] - V[v1]).T)
            rhs[v0] -= b.T
            rhs[v1] += b.T

    return rhs


def gauss_thinning(
    V,
    F,
    num_iterations=100,
    min_cone_angle=2.5,
    smooth=1e-5,
    start_angle=25,
    radius=0.1,
    sigma=2,
    render=False,
):
    cone_angle = start_angle
    eps = 1e-3

    nv = len(V)
    center(V)

    TT = triangle_adjacency(F, nv)
    C = igl.cotmatrix_entries(V, F)
    L = igl.cotmatrix(V, F)
    M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_BARYCENTRIC).todense()

    if smooth:
        A = -L + smooth * np.dot(L.T, L) + eps * M
    else:
        A = -L + eps * M

    V_dash = np.copy(V)

    if render:
        plot = mp.Viewer({})
        mesh = plot.add_mesh(V_dash, F, c=V_dash[:, 0])
        display(plot._renderer)  # type: ignore

    for i in range(num_iterations):
        N = igl.per_face_normals(V_dash, F, np.array([1, 0, 0], dtype=float))
        B = igl.barycenter(V_dash, F)

        nbhs = collectNeighbors(TT, B, N, radius, cone_angle)
        if cone_angle > min_cone_angle:
            cone_angle *= 0.95

        N2 = fit_normals(nbhs, N, cone_angle, sigma)
        rot = find_rotations(N, N2)
        b = assemble_RHS(C, V_dash, F, rot)

        V_dash = np.asarray(LA.solve(A, eps * np.dot(M, V_dash) - b))

        if render:
            plot.update_object(oid=mesh, vertices=V_dash, colors=V_dash[:, 0], faces=F)

        print(f"Finished iteration {i}.")

    return V_dash
