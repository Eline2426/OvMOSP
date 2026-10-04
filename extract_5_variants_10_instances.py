import argparse
import math
import os
import random
from collections import deque

BOROUGHS = {
    "Manhattan":    (40.7580, -73.9855),
    "Brooklyn":     (40.6782, -73.9442),
    "Queens":       (40.7282, -73.7949),
    "Bronx":        (40.8448, -73.8648),
    "StatenIsland": (40.5795, -74.1502),
}


def read_gr(fname):
    n = m = None
    arcs = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("c"):
                continue
            parts = line.split()
            if parts[0] == "p":
                n = int(parts[2])
                m = int(parts[3])
            elif parts[0] == "a":
                arcs.append((int(parts[1]), int(parts[2]), int(parts[3])))
    return n, arcs


def read_coords(fname):
    coords = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("c") or line.startswith("p"):
                continue
            parts = line.split()
            if parts[0] == "v":
                node_id = int(parts[1])
                lon = int(parts[2]) / 1e6
                lat = int(parts[3]) / 1e6
                coords[node_id] = (lat, lon)
    return coords


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlmb = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlmb / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


def merge_all(distance_file, time_file, extra_files):
    n_d, arcs_d = read_gr(distance_file)
    n_t, arcs_t = read_gr(time_file)
    extra = [read_gr(f)[1] for f in extra_files]

    all_lists = [arcs_d, arcs_t] + extra

    if len({len(a) for a in all_lists}) != 1:
        raise ValueError("The .gr files do not have the same number of arcs")

    merged = []
    for rows in zip(*all_lists):
        u0, v0, _ = rows[0]
        if any((r[0], r[1]) != (u0, v0) for r in rows):
            raise ValueError("The arcs are not aligned across the .gr files")
        costs = [w for (_, _, w) in rows]
        merged.append((u0, v0, costs[0], tuple(costs[1:])))
    return n_d, merged


def largest_component(nodes, arcs):
    node_set = set(nodes)
    adj = {}
    for u, v, d, c in arcs:
        if u in node_set and v in node_set:
            adj.setdefault(u, set()).add(v)
            adj.setdefault(v, set()).add(u)

    seen_global = set()
    best_component = set()
    for start in node_set:
        if start in seen_global or start not in adj:
            continue
        seen = {start}
        q = deque([start])
        while q:
            x = q.popleft()
            for y in adj.get(x, []):
                if y not in seen:
                    seen.add(y)
                    q.append(y)
        seen_global |= seen
        if len(seen) > len(best_component):
            best_component = seen

    return best_component


def write_instance(fname, n, arcs, s, t, p):
    with open(fname, "w") as f:
        f.write(f"{n} {len(arcs)} {p} {s} {t}\n")
        for u, v, d, c in arcs:
            f.write(f"{u} {v} {d} " + " ".join(map(str, c)) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--distance", required=True)
    ap.add_argument("--time", required=True)
    ap.add_argument("--elevation", required=True)
    ap.add_argument("--avg_degree", default=None)
    ap.add_argument("--hop_count", default=None)
    ap.add_argument("--coords", required=True)
    ap.add_argument("--radius_km", type=float, default=8.0)
    ap.add_argument("--nb_instances", type=int, default=10)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    extra_files = [args.elevation]
    if args.avg_degree:
        extra_files.append(args.avg_degree)
    if args.hop_count:
        extra_files.append(args.hop_count)

    print("Fusion des fichiers (une seule fois pour tout le reseau)...")
    n, arcs = merge_all(args.distance, args.time, extra_files)
    p = len(extra_files) + 1
    print(f"Reseau complet : n={n}, m={len(arcs)}, p={p}")

    print("Lecture des coordonnees...")
    coords = read_coords(args.coords)

    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)

    summary = []

    for borough, (lat0, lon0) in BOROUGHS.items():
        print(f"\n=== {borough} (centre {lat0:.4f}, {lon0:.4f}, rayon {args.radius_km} km) ===")

        nodes_in_radius = {
            nid for nid, (lat, lon) in coords.items()
            if haversine_km(lat0, lon0, lat, lon) <= args.radius_km
        }

        component = largest_component(nodes_in_radius, arcs)
        if len(component) < 10:
            print(f"  [avertissement] composante trop petite ({len(component)} noeuds) -- "
                  f"augmentez --radius_km pour {borough}.")
            continue

        sub_arcs = [(u, v, d, c) for u, v, d, c in arcs if u in component and v in component]
        component_list = sorted(component)
        remap = {old: i for i, old in enumerate(component_list)}
        new_arcs = [(remap[u], remap[v], d, c) for u, v, d, c in sub_arcs]

        print(f"  Sous-reseau : n={len(component)}, m={len(new_arcs)}")

        adj = {}
        for u, v, d, c in new_arcs:
            adj.setdefault(u, []).append(v)

        def bfs_dist(seed_node):
            dist = {seed_node: 0}
            q = deque([seed_node])
            while q:
                x = q.popleft()
                for y in adj.get(x, []):
                    if y not in dist:
                        dist[y] = dist[x] + 1
                        q.append(y)
            return dist

        nb_generated = 0
        attempts = 0
        min_hops = 10

        while nb_generated < args.nb_instances and attempts < args.nb_instances * 50:
            attempts += 1

            if attempts % (args.nb_instances * 10) == 0 and min_hops > 2:
                min_hops -= 2
                print(f"    [ajustement] min_hops reduit a {min_hops} "
                      f"(peu de candidats trouves)")

            s_node = rng.choice(component_list)
            s_new = remap[s_node]
            dist_map = bfs_dist(s_new)
            candidates = [v for v, dd in dist_map.items() if dd >= min_hops]
            if not candidates:
                continue
            t_new = rng.choice(candidates)

            nb_generated += 1
            out_fname = os.path.join(args.out_dir, f"{borough}_q{nb_generated:02d}.txt")
            write_instance(out_fname, len(component), new_arcs, s_new, t_new, p)
            print(f"  [{nb_generated}/{args.nb_instances}] s={s_new} t={t_new} "
                  f"(dist topo={dist_map[t_new]}) -> {out_fname}")
            summary.append((borough, nb_generated, len(component), len(new_arcs), s_new, t_new))

        if nb_generated < args.nb_instances:
            print(f"  [avertissement] seulement {nb_generated}/{args.nb_instances} "
                  f"paires generees pour {borough} -- reseau peut-etre trop petit.")

    print(f"\n{'=' * 60}")
    print(f"Termine : {len(summary)} instances generees au total dans {args.out_dir}/")
    print(f"{'=' * 60}")
