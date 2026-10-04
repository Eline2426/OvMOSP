import time
import os
import sys
import random
from heapq import heappush, heappop

import cplex

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False


N_SCALARIZATIONS = 50


RANDOM_SEED = 12345


def dominates(a, b):
    return all(x <= y for x, y in zip(a, b)) and any(
        x < y for x, y in zip(a, b)
    )


def update_UBC(UBC, UBC_set, new_vec):
    if new_vec in UBC_set:
        return UBC, UBC_set

    for v in UBC:
        if dominates(v, new_vec):
            return UBC, UBC_set

    new_UBC = [v for v in UBC if not dominates(new_vec, v)]
    new_UBC.append(new_vec)
    new_UBC_set = set(new_UBC)

    return new_UBC, new_UBC_set


def build_arc_index(arcs):
    index = {}
    for (u, v, d, c) in arcs:
        if (u, v) not in index:
            index[(u, v)] = (d, c)
    return index


def path_costs(path, arc_index, p):
    d = 0
    c = [0] * p

    for i in range(len(path) - 1):
        dd, cc = arc_index[(path[i], path[i + 1])]
        d += dd
        c = [x + y for x, y in zip(c, cc)]

    return d, tuple(c)


def dijkstra_forward(adj_out, n, s):
    dist = [float('inf')] * n
    pred = [-1] * n
    dist[s] = 0

    heap = [(0.0, s)]

    while heap:
        du, u = heappop(heap)

        if du > dist[u]:
            continue

        for (v, w, _cc) in adj_out[u]:
            nd = du + w
            if nd < dist[v]:
                dist[v] = nd
                pred[v] = u
                heappush(heap, (nd, v))

    return dist, pred


def reconstruct_from_pred(pred, s, t):
    if pred[t] == -1 and s != t:
        return None

    path = []
    cur = t

    while cur != -1:
        path.append(cur)
        if cur == s:
            break
        cur = pred[cur]

    if not path or path[-1] != s:
        return None

    path.reverse()

    return path


def dijkstra_backward_dist(adj_in, n, t):
    dist = [float('inf')] * n
    dist[t] = 0.0

    heap = [(0.0, t)]

    while heap:
        dv, v = heappop(heap)

        if dv > dist[v]:
            continue

        for (u, w, _cc) in adj_in[v]:
            nd = dv + w
            if nd < dist[u]:
                dist[u] = nd
                heappush(heap, (nd, u))

    return dist


def dijkstra_backward_vector(adj_in, n, t, p):
    LB = [[float('inf')] * p for _ in range(n)]

    for k in range(p):
        dist_k = [float('inf')] * n
        dist_k[t] = 0.0
        heap = [(0.0, t)]

        while heap:
            dv, v = heappop(heap)

            if dv > dist_k[v]:
                continue

            for (u, _d, cc) in adj_in[v]:
                nd = dv + cc[k]
                if nd < dist_k[u]:
                    dist_k[u] = nd
                    heappush(heap, (nd, u))

        for v in range(n):
            LB[v][k] = dist_k[v]

    return LB


_bes_cache = {}


def test_BES(arcs, n, s, t, p, path_P, arc_index, timelimit=60.0):
    _, critP = path_costs(path_P, arc_index, p)

    key = (id(arcs), critP, n, s, t, p)
    if key in _bes_cache:
        return _bes_cache[key]

    nb = len(arcs)

    model = cplex.Cplex()
    model.set_results_stream(None)
    model.parameters.timelimit.set(timelimit)

    x_names = [f"x{i}" for i in range(nb)]
    v_names = [f"v{i}" for i in range(p)]

    model.variables.add(
        names=x_names,
        types=["B"] * nb,
        lb=[0] * nb,
        ub=[1] * nb
    )

    model.variables.add(
        names=v_names,
        lb=[0] * p,
        ub=[cplex.infinity] * p
    )

    obj = [(v_names[i], 1.0) for i in range(p)]
    model.objective.set_sense(model.objective.sense.maximize)
    model.objective.set_linear(obj)


    out = [[] for _ in range(n)]
    inc = [[] for _ in range(n)]

    for i, (u, v, _, _) in enumerate(arcs):
        out[u].append(i)
        inc[v].append(i)

    for node in range(n):
        ind = [x_names[i] for i in out[node]] + [x_names[i] for i in inc[node]]
        val = [1] * len(out[node]) + [-1] * len(inc[node])

        if not ind:
            continue

        if node == s:
            rhs = 1
        elif node == t:
            rhs = -1
        else:
            rhs = 0

        model.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=ind, val=val)],
            senses=["E"],
            rhs=[rhs]
        )


    for k in range(p):
        coeff = [arcs[i][3][k] for i in range(nb)]
        ind = x_names + [v_names[k]]
        val = coeff + [1.0]

        model.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=ind, val=val)],
            senses=["L"],
            rhs=[critP[k]]
        )

    try:
        model.solve()

        sample_mem()


        if model.solution.get_status() not in (1, 101, 102):
            return {"status": "error", "path": None, "v": None, "sum_v": None}

        if not model.solution.is_primal_feasible():
            result = {"status": "infeasible", "path": None, "v": None, "sum_v": None}
            _bes_cache[key] = result
            return result

        sol_x = [model.solution.get_values(xi) for xi in x_names]
        sol_v = [model.solution.get_values(vi) for vi in v_names]

        sum_v = sum(sol_v)
        tolerance = 1e-7

        if sum_v <= tolerance:
            result = {"status": "efficient", "path": path_P[:], "v": sol_v, "sum_v": sum_v}
        else:
            dist_Q = sum(arcs[i][2] * sol_x[i] for i in range(nb))
            result = {"status": "dominated", "sol_x": sol_x, "v": sol_v,
                      "sum_v": sum_v, "dist": dist_Q}

        _bes_cache[key] = result
        return result

    except Exception:
        return {"status": "error", "path": None, "v": None, "sum_v": None}


def reconstruct_path(arcs, sol_x, s, t, th=0.5):
    succ = {}

    for (u, v, _, _), val in zip(arcs, sol_x):
        if val > th:
            succ[u] = v

    if s not in succ:
        return None

    path = [s]
    cur = s
    visited = set()

    while cur != t:
        if cur in visited:
            return None

        visited.add(cur)

        if cur not in succ:
            return None

        cur = succ[cur]
        path.append(cur)

    return path


def generate_weight_vectors(p, n_scalarizations=50, seed=12345):
    if p <= 0:
        return []

    if n_scalarizations <= 0:
        return []

    rng = random.Random(seed)

    weights = []

    for _ in range(n_scalarizations):
        values = [rng.expovariate(1.0) for _ in range(p)]
        total = sum(values)
        lam = tuple(x / total for x in values)
        weights.append(lam)

    return weights


def scalar_path(adj_out, n, s, t, weights):
    p = len(weights)
    dist = [float('inf')] * n
    pred = [-1] * n
    dist[s] = 0.0

    heap = [(0.0, s)]

    while heap:
        du, u = heappop(heap)

        if du > dist[u]:
            continue

        if u == t:
            break

        for (v, _d, cc) in adj_out[u]:
            w = sum(weights[k] * cc[k] for k in range(p))
            nd = du + w

            if nd < dist[v]:
                dist[v] = nd
                pred[v] = u
                heappush(heap, (nd, v))

    if dist[t] == float('inf'):
        return None

    return reconstruct_from_pred(pred, s, t)


def scalar_path_with_cost(adj_out, n, s, t, weights):
    p = len(weights)
    dist = [float('inf')] * n
    pred = [-1] * n
    cost_vec = [None] * n
    dist[s] = 0.0
    cost_vec[s] = tuple([0.0] * p)

    heap = [(0.0, s)]

    while heap:
        du, u = heappop(heap)

        if du > dist[u]:
            continue

        if u == t:
            break

        for (v, _d, cc) in adj_out[u]:
            w = sum(weights[k] * cc[k] for k in range(p))
            nd = du + w

            if nd < dist[v]:
                dist[v] = nd
                pred[v] = u
                cost_vec[v] = tuple(x + y for x, y in zip(cost_vec[u], cc))
                heappush(heap, (nd, v))

    if dist[t] == float('inf'):
        return None, None

    path = reconstruct_from_pred(pred, s, t)
    return path, cost_vec[t]


def compute_normalization_scales(adj_out, n, s, t, p):
    scales = [1.0] * p
    for k in range(p):
        w = [0.0] * p
        w[k] = 1.0
        path, cvec = scalar_path_with_cost(adj_out, n, s, t, w)
        if path and cvec is not None and cvec[k] > 0:
            scales[k] = cvec[k]
    return scales


_PEAK = [0.0]


def get_rss_mb():
    if not _HAS_PSUTIL:
        return None
    return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)


def sample_mem():
    rss = get_rss_mb()
    if rss is not None and rss > _PEAK[0]:
        _PEAK[0] = rss
    return rss


def get_peak_os_mb():
    if _HAS_PSUTIL:
        mi = psutil.Process(os.getpid()).memory_info()
        if hasattr(mi, "peak_wset"):
            return mi.peak_wset / (1024 * 1024)

    try:
        import resource
        r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == "darwin":
            return r / (1024 * 1024)
        return r / 1024
    except Exception:
        return None


def _fmt(x):
    return "NA" if x is None else f"{x:.1f}"


def structures_size_estimate(stock, UBC, stack_len):
    total_stock = sum(len(s) for s in stock)
    return {
        "stock_vectors_total": total_stock,
        "UBC_vectors": len(UBC),
        "stack_depth": stack_len
    }


def initialization(arcs, n, s, t, p, arc_index, adj_out,
                   scalarization_count=50, random_seed=12345,
                   bet_timelimit=60.0):


    UBC = []
    UBC_set = set()
    best_dist = {}
    solutions = []
    UB = float('inf')

    nb_BES_calls = 0

    t_init_start = time.perf_counter()

    scalarization_time = 0.0
    UB_before_scalarization = float('inf')
    UB_after_scalarization = float('inf')
    bound_improvement_scalarization = 0.0

    ddist, pred = dijkstra_forward(adj_out, n, s)

    if ddist[t] == float('inf'):
        init_time = time.perf_counter() - t_init_start
        return (UBC, UBC_set, UB, solutions, best_dist, 0, False, None,
                nb_BES_calls, init_time, scalarization_time,
                UB_before_scalarization, UB_after_scalarization,
                bound_improvement_scalarization, 0)

    P0 = reconstruct_from_pred(pred, s, t)
    d0, c0 = path_costs(P0, arc_index, p)

    print("\n========== ETAPE 1 : PLUS COURT CHEMIN SELON d ==========")
    print(f"P0 = {P0}  d(P0) = {d0}  c(P0) = {c0}")

    print("\n========== ETAPE 2 : TEST D'EFFICACITE DE BENSON (BET) ==========")

    bes = test_BES(arcs, n, s, t, p, P0, arc_index, timelimit=bet_timelimit)
    nb_BES_calls += 1

    if bes["status"] == "efficient":
        print("P0 est efficace -> arret immediat.")

        UBC.append(c0)
        UBC_set.add(c0)
        best_dist[c0] = d0
        solutions.append(P0[:])

        init_time = time.perf_counter() - t_init_start

        return (UBC, UBC_set, d0, solutions, best_dist, 1, True, P0,
                nb_BES_calls, init_time, scalarization_time,
                d0, d0, 0.0, 0)

    elif bes["status"] == "dominated":
        print("P0 n'est pas efficace ; recuperation du chemin dominant.")

        q = reconstruct_path(arcs, bes["sol_x"], s, t)
        if q is None:
            raise RuntimeError("Impossible de reconstruire le chemin fourni par BET.")

        dq, cq = path_costs(q, arc_index, p)
        print(f"x = {q}  d(x) = {dq}  c(x) = {cq}")

        UB = dq
        UBC, UBC_set = update_UBC(UBC, UBC_set, cq)
        best_dist[cq] = dq
        solutions.append(q[:])

        print(f"Initialisation : UB = {UB}  UBC = {UBC}")

    else:
        raise RuntimeError("Le test BET n'a pas pu etre effectue correctement.")

    UB_before_scalarization = UB

    print("\n========== ETAPE 3 : ENRICHISSEMENT (PHASE 2 : SCALARISATION ALEATOIRE) ==========")

    if scalarization_count <= 0:
        print("Scalarisation DESACTIVEE -> alpha = 0, Phase 2 sautee.")
        weights = []
        scalarization_time = 0.0

    else:
        scales = compute_normalization_scales(adj_out, n, s, t, p)
        print(f"Echelles de normalisation par critere : {scales}")

        weights = generate_weight_vectors(
            p,
            n_scalarizations=scalarization_count,
            seed=random_seed
        )

        print(f"Nombre de vecteurs lambda aleatoires : {len(weights)}")
        print(f"Graine aleatoire utilisee : {random_seed}")

        if weights:
            print("Exemple des premiers vecteurs lambda :")
            for i in range(min(5, len(weights))):
                print(f"  lambda {i + 1} : {weights[i]}  somme = {sum(weights[i]):.12f}")

        t_scalar_start = time.perf_counter()

        for i, w in enumerate(weights, start=1):
            w_normalized = tuple(w[k] / scales[k] for k in range(p))

            path = scalar_path(adj_out, n, s, t, w_normalized)

            if not path:
                continue

            d, c = path_costs(path, arc_index, p)

            old_len = len(UBC)
            UBC, UBC_set = update_UBC(UBC, UBC_set, c)

            if c not in best_dist or d < best_dist[c]:
                best_dist[c] = d
                solutions.append(path[:])

            if d < UB:
                UB = d

            if len(UBC) > old_len:
                print(f"  Scalarisation {i}/{len(weights)} -> nouveau vecteur efficace {c}")

        scalarization_time = time.perf_counter() - t_scalar_start


    alpha = len(weights)

    UB_after_scalarization = UB
    bound_improvement_scalarization = UB_before_scalarization - UB_after_scalarization

    nb_initial = len(UBC)

    init_time = time.perf_counter() - t_init_start

    print(f"\nFin de l'initialisation. UB = {UB}  |UBC| = {len(UBC)}")
    print(f"Alpha (nb de vecteurs de scalarisation utilises) : {alpha}")
    print(f"Temps de la phase de scalarisation : {scalarization_time:.4f} s")
    print(f"UB avant scalarisation : {UB_before_scalarization}")
    print(f"UB apres scalarisation : {UB_after_scalarization}")
    print(f"Borne amelioree grace a la scalarisation (delta UB) : {bound_improvement_scalarization}")

    return (UBC, UBC_set, UB, solutions, best_dist, nb_initial, False, P0,
            nb_BES_calls, init_time, scalarization_time,
            UB_before_scalarization, UB_after_scalarization,
            bound_improvement_scalarization, alpha)


def branch_and_bound_fast(arcs, n, s, t, p,
                          scalarization_count=50,
                          random_seed=12345,
                          bet_timelimit=60.0,
                          track_memory=True,
                          memory_sample_every=2000):


    start_time = time.perf_counter()

    _PEAK[0] = 0.0
    baseline_mb = sample_mem() if track_memory else None

    arc_index = build_arc_index(arcs)

    adj_out = [[] for _ in range(n)]
    adj_in = [[] for _ in range(n)]

    for (u, v, d, cc) in arcs:
        adj_out[u].append((v, d, cc))
        adj_in[v].append((u, d, cc))


    for u in range(n):
        adj_out[u].sort(key=lambda x: x[1])

    (UBC, UBC_set, UB, solutions, best_dist, nb_initial, stop, P0,
     nb_BES_init, init_time, scalarization_time,
     UB_before_scalarization, UB_after_scalarization,
     bound_improvement_scalarization, alpha) = initialization(
        arcs, n, s, t, p, arc_index, adj_out,
        scalarization_count=scalarization_count,
        random_seed=random_seed,
        bet_timelimit=bet_timelimit
    )

    if track_memory:
        sample_mem()

    if UB == float('inf'):
        return None


    if stop:
        elapsed = time.perf_counter() - start_time
        d0, c0 = path_costs(P0, arc_index, p)

        if track_memory:
            sample_mem()

        return {
            "best_path": P0, "best_d": d0, "best_c": c0, "UBC": UBC,
            "cpu": elapsed, "nodes": 0, "pruned": 0,
            "pruned_rule1": 0, "pruned_rule2": 0, "pruned_rule3": 0,
            "calls_BES": nb_BES_init,
            "nb_initial": 1, "nb_added": 0, "total_efficient": 1,
            "early_stop": True,
            "peak_memory_mb": _PEAK[0] if track_memory else None,
            "peak_memory_os_mb": get_peak_os_mb(),
            "baseline_memory_mb": baseline_mb,
            "max_stock_vectors": 0,
            "max_UBC": 1,
            "bes_cache_entries": len(_bes_cache),
            "structures_peak": {"stock_vectors_total": 0, "UBC_vectors": 1, "stack_depth": 0},
            "nb_labels_generated": 0,
            "nb_labels_stored": 0,
            "init_time": init_time,
            "scalarization_time": scalarization_time,
            "UB_before_scalarization": UB_before_scalarization,
            "UB_after_scalarization": UB_after_scalarization,
            "bound_improvement_scalarization": bound_improvement_scalarization,
            "tree_depth": 0,
            "alpha": alpha,
            "random_seed": random_seed,
            "n_scalarizations": scalarization_count,
        }

    LB_dist = dijkstra_backward_dist(adj_in, n, t)
    LB_vec = dijkstra_backward_vector(adj_in, n, t, p)


    local_dom = [[] for _ in range(n)]

    nb_nodes = 0
    nb_pruned_rule1 = 0
    nb_pruned_rule2 = 0
    nb_pruned_rule3 = 0
    nb_BES = 0
    nb_efficient_added = 0
    max_stack_depth = 0

    stock_total = 0
    max_stock_total = 0
    max_UBC = len(UBC)

    def dfs(node, path, cost_vec, dval):
        nonlocal UB, UBC, UBC_set, best_dist, solutions
        nonlocal nb_nodes, nb_pruned_rule1, nb_pruned_rule2, nb_pruned_rule3
        nonlocal nb_BES, nb_efficient_added, max_stack_depth
        nonlocal stock_total, max_stock_total, max_UBC

        nb_nodes += 1
        path.append(node)

        max_stack_depth = max(max_stack_depth, len(path))

        if track_memory and nb_nodes % memory_sample_every == 0:
            sample_mem()


        for v in local_dom[node]:
            if dominates(v, cost_vec):
                path.pop()
                nb_pruned_rule2 += 1
                return

        new_local = [cost_vec]
        for v in local_dom[node]:
            if not dominates(cost_vec, v):
                new_local.append(v)
        stock_total += len(new_local) - len(local_dom[node])

        local_dom[node] = new_local

        if stock_total > max_stock_total:
            max_stock_total = stock_total


        future_min = tuple(cost_vec[k] + LB_vec[node][k] for k in range(p))

        if any(dominates(u, future_min) for u in UBC):
            path.pop()
            nb_pruned_rule3 += 1
            return


        if node == t:

            dominated_by_UBC = any(dominates(u, cost_vec) for u in UBC)

            if not dominated_by_UBC:
                nb_BES += 1

                bes = test_BES(arcs, n, s, t, p, path, arc_index,
                               timelimit=bet_timelimit)

                if bes["status"] == "dominated":
                    q = reconstruct_path(arcs, bes["sol_x"], s, t)

                    if q:
                        dq, cq = path_costs(q, arc_index, p)

                        old_len = len(UBC)
                        UBC, UBC_set = update_UBC(UBC, UBC_set, cq)

                        if len(UBC) > old_len:
                            nb_efficient_added += 1

                        if cq not in best_dist or dq < best_dist[cq]:
                            best_dist[cq] = dq
                            solutions.append(q[:])

                        if dq < UB:
                            UB = dq

                elif bes["status"] == "efficient":
                    d_curr = dval
                    c_curr = tuple(cost_vec)

                    old_len = len(UBC)
                    UBC, UBC_set = update_UBC(UBC, UBC_set, c_curr)

                    if len(UBC) > old_len:
                        nb_efficient_added += 1

                    if c_curr not in best_dist or d_curr < best_dist[c_curr]:
                        best_dist[c_curr] = d_curr
                        solutions.append(path[:])

                    if d_curr < UB:
                        UB = d_curr

                else:


                    raise RuntimeError("BET not solved to optimality (time limit?)")

            if len(UBC) > max_UBC:
                max_UBC = len(UBC)

            path.pop()
            return


        for v, d_uv, c_uv in adj_out[node]:

            if v in path:
                continue

            nd = dval + d_uv


            if nd + LB_dist[v] > UB:
                nb_pruned_rule1 += 1
                continue

            ncost = tuple(x + y for x, y in zip(cost_vec, c_uv))

            dfs(v, path, ncost, nd)

        path.pop()

    dfs(s, [], tuple([0] * p), 0)


    best_path = None
    best_d = float('inf')
    best_c = None

    for path in solutions:
        d, c = path_costs(path, arc_index, p)

        if c in UBC_set:
            if d < best_d:
                best_d = d
                best_path = path[:]
                best_c = c

    if best_path is None and best_dist:
        candidates = [(c, d) for c, d in best_dist.items() if c in UBC_set]

        if candidates:
            best_c, best_d = min(candidates, key=lambda x: x[1])

            for path in solutions:
                _, c = path_costs(path, arc_index, p)
                if c == best_c:
                    best_path = path[:]
                    break

    elapsed = time.perf_counter() - start_time
    total_efficient = len(UBC)

    nb_pruned = nb_pruned_rule1 + nb_pruned_rule2 + nb_pruned_rule3

    nb_labels_stored = sum(len(s) for s in local_dom)

    if track_memory:
        sample_mem()

    structures_peak = structures_size_estimate(local_dom, UBC, max_stack_depth)

    return {
        "best_path": best_path,
        "best_d": best_d,
        "best_c": best_c,
        "UBC": UBC,
        "cpu": elapsed,
        "nodes": nb_nodes,
        "pruned": nb_pruned,
        "pruned_rule1": nb_pruned_rule1,
        "pruned_rule2": nb_pruned_rule2,
        "pruned_rule3": nb_pruned_rule3,
        "calls_BES": nb_BES + nb_BES_init,
        "nb_initial": nb_initial,
        "nb_added": nb_efficient_added,
        "total_efficient": total_efficient,
        "early_stop": False,
        "peak_memory_mb": _PEAK[0] if track_memory else None,
        "peak_memory_os_mb": get_peak_os_mb(),
        "baseline_memory_mb": baseline_mb,
        "max_stock_vectors": max_stock_total,
        "max_UBC": max_UBC,
        "bes_cache_entries": len(_bes_cache),
        "structures_peak": structures_peak,
        "nb_labels_generated": nb_nodes,
        "nb_labels_stored": nb_labels_stored,
        "init_time": init_time,
        "scalarization_time": scalarization_time,
        "UB_before_scalarization": UB_before_scalarization,
        "UB_after_scalarization": UB_after_scalarization,
        "bound_improvement_scalarization": bound_improvement_scalarization,
        "tree_depth": max_stack_depth,
        "alpha": alpha,
        "random_seed": random_seed,
        "n_scalarizations": scalarization_count,
    }


def read_instance(fname):
    with open(fname, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    first = lines[0].split()

    n = int(first[0])
    m = int(first[1])
    p = int(first[2])

    s = int(first[3]) if len(first) > 3 else 0
    t = int(first[4]) if len(first) > 4 else n - 1

    arcs = []
    min_node = float("inf")
    max_node = -float("inf")

    for line in lines[1:1 + m]:
        values = line.split()

        u = int(values[0])
        v = int(values[1])
        d = float(values[2])
        c = tuple(float(x) for x in values[3:])

        arcs.append((u, v, d, c))

        min_node = min(min_node, u, v)
        max_node = max(max_node, u, v)


    if min_node == 1:
        arcs = [(u - 1, v - 1, d, c) for u, v, d, c in arcs]
        s -= 1
        t -= 1
        n = int(max_node)
    else:
        n = int(max_node + 1)

    return arcs, n, len(arcs), p, s, t


if __name__ == "__main__":

    fname = sys.argv[1] if len(sys.argv) > 1 else "instance.txt"

    arcs, n, m, p, s, t = read_instance(fname)

    print(f"\nInstance : {fname}\nNoeuds = {n}\nArcs = {m}"
          f"\nCriteres secondaires = {p}\nSource = {s}\nCible = {t}\n")

    if not _HAS_PSUTIL:
        print("ATTENTION : psutil n'est pas installe, la memoire "
              "RSS ne sera pas mesuree (pip install psutil).\n")


    if len(sys.argv) > 2:
        arg = sys.argv[2].strip().lower()
        if arg in ("0", "off", "none", "disable", "disabled"):
            scalarization_count = 0
        else:
            scalarization_count = int(arg)
    else:
        scalarization_count = N_SCALARIZATIONS


    random_seed = int(sys.argv[3]) if len(sys.argv) > 3 else RANDOM_SEED


    bet_timelimit = float(sys.argv[4]) if len(sys.argv) > 4 else 6000.0

    print(f"Nombre de scalarisations aleatoires = {scalarization_count}")
    print(f"Seed aleatoire = {random_seed}")
    print("\nLes vecteurs lambda satisfont :")
    print("lambda_i > 0 et somme(lambda_i) = 1")

    res = branch_and_bound_fast(
        arcs, n, s, t, p,
        scalarization_count=scalarization_count,
        random_seed=random_seed,
        bet_timelimit=bet_timelimit,
        track_memory=True,
        memory_sample_every=2000
    )

    if res:
        print("\n================================================")
        print("RESULTATS FINAUX")
        print("================================================")
        print(f"Nombre de scalarisations aleatoires : {res['n_scalarizations']}")
        print(f"Alpha (nb vecteurs de scalarisation) : {res['alpha']}")
        print(f"Seed utilise : {res['random_seed']}")
        print(f"Noeuds explores (labels generes) : {res['nodes']}")
        print(f"  dont stockes dans stock[v] : {res['nb_labels_stored']}")
        print(f"Profondeur de l'arborescence : {res['tree_depth']}")
        print(f"Elagages totaux : {res['pruned']}")
        print(f"  Regle 1 (borne scalaire) : {res['pruned_rule1']}")
        print(f"  Regle 2 (dominance locale) : {res['pruned_rule2']}")
        print(f"  Regle 3 (borne vectorielle) : {res['pruned_rule3']}")
        print(f"Appels a BET (total) : {res['calls_BES']}")
        print(f"Temps CPU (total) : {res['cpu']:.2f} s")
        print(f"  dont temps d'initialisation : {res['init_time']:.4f} s")
        print(f"  dont temps de scalarisation : {res['scalarization_time']:.4f} s")
        print(f"UB avant scalarisation : {res['UB_before_scalarization']}")
        print(f"UB apres scalarisation : {res['UB_after_scalarization']}")
        print(f"Borne amelioree par scalarisation (delta UB) : {res['bound_improvement_scalarization']}")
        print(f"Memoire de base (avant resolution) : {_fmt(res['baseline_memory_mb'])} Mo")
        print(f"Pic de memoire echantillonne (RSS)  : {_fmt(res['peak_memory_mb'])} Mo")
        print(f"Pic de memoire (systeme)            : {_fmt(res['peak_memory_os_mb'])} Mo")
        print(f"Taille max de sum|stock[v]|         : {res['max_stock_vectors']}")
        print(f"|UBC| max                           : {res['max_UBC']}")
        print(f"Entrees du cache BET                : {res['bes_cache_entries']}")
        print(f"Structures (pic) : {res['structures_peak']}")
        print(f"Solutions efficaces initiales : {res['nb_initial']}")
        print(f"Solutions efficaces ajoutees : {res['nb_added']}")
        print(f"Nombre total de solutions efficaces : {res['total_efficient']}")
        print(f"Meilleur vecteur de criteres secondaires : {res['best_c']}")
        print(f"Distance minimale parmi les efficaces : {res['best_d']}")
        print(f"Meilleur chemin : {res['best_path']}")
        print(f"Arret immediat apres BET (phase 1) : {res['early_stop']}")


        print("RESULT_LINE: " + ";".join([
            f"alpha={res['alpha']}",
            f"UB_scalarisation={res['UB_after_scalarization']}",
            f"UB_final={res['best_d']}",
            f"cpu_scalarisation={res['scalarization_time']:.4f}",
            f"cpu_total={res['cpu']:.4f}",
            f"profondeur_arbre={res['tree_depth']}",
            f"peak_mem_mb={_fmt(res['peak_memory_mb'])}",
            f"peak_mem_os_mb={_fmt(res['peak_memory_os_mb'])}",
            f"base_mem_mb={_fmt(res['baseline_memory_mb'])}",
            f"max_stock={res['max_stock_vectors']}",
            f"max_UBC={res['max_UBC']}",
            f"bet_cache={res['bes_cache_entries']}",
        ]))

    else:
        print("Aucun chemin trouve.")
        print("RESULT_LINE: status=NO_PATH")
