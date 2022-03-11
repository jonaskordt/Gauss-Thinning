import math
import numpy as np
from numpy import linalg as LA
import igl
from gauss_thinning import *

# add custom pybind lib
import sys
sys.path.insert(0,'../third_party/DevelopableApproximationViaGaussImageThinning')
from mainParallel import getMassmatrix, triangleAdjacency, createInitialNBHS, createMatrixXd, createRotations, getCotmatrixEntries, thinningIteration, getCotmatrix


# TODO: Implement this in C++
def expand_path(path, r, V, F):
    nv = len(V)
    adj = triangle_adjacency(F, nv)
    B = igl.barycenter(V, F)
    stack = []
    flag = [-1 for _ in range(len(B))]
    flag_result = [False for _ in range(len(B))]

    for i in path:
        stack.append(i)
        flag[i] = i

        while len(stack) > 0:
            id = stack.pop()

            
            flag_result[id] = True

            for j in adj[id]:
                if flag[j] != i and LA.norm(B[i] - B[j]) < r:
                    stack.append(j)
                    flag[j] = i

    return flag_result


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
    r_squared = radius**2

    nv = len(V)

    M = getMassmatrix(V, F)
    L = getCotmatrix(V, F)
    N = createMatrixXd()
    N2 = createMatrixXd()
    B = createMatrixXd()
    tt = triangleAdjacency(F, nv)
    nbhs = createInitialNBHS(F)
    rot = createRotations()
    C = getCotmatrixEntries(V, F)
    b = createMatrixXd()

    V_dash = np.copy(V)
    
    for i in range(num_iterations):
        V_dash = thinningIteration(
            V_dash,
            F,
            active_triangles,
            M,
            L,
            N,
            N2,
            B,
            tt,
            nbhs,
            rot,
            C,
            b,
            r_squared,
            cone_angle,
            min_cone_angle,
            sigma,
            eps,
            smooth
        )

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