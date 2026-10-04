import argparse
import glob
import os


def truncate_instance(input_path, output_path, remove_last_n):
    with open(input_path, "r") as f:
        lines = [line for line in f if line.strip()]

    header_tokens = lines[0].split()
    n, m, p, s, t = header_tokens
    p = int(p)
    new_p = p - remove_last_n

    if new_p < 0:
        raise ValueError(f"{input_path} : impossible de retirer {remove_last_n} colonnes, p={p}")

    new_lines = [f"{n} {m} {new_p} {s} {t}\n"]

    for line in lines[1:]:
        tokens = line.split()
        u, v, d = tokens[0], tokens[1], tokens[2]
        c = tokens[3:3 + p]
        if len(c) != p:
            raise ValueError(f"{input_path} : ligne avec {len(c)} couts, attendu p={p} -> {line.strip()}")
        c_kept = c[: len(c) - remove_last_n]
        new_lines.append(" ".join([u, v, d] + c_kept) + "\n")

    with open(output_path, "w") as f:
        f.writelines(new_lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True)
    ap.add_argument("--out_dir_p3", required=True)
    ap.add_argument("--out_dir_p2", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir_p3, exist_ok=True)
    os.makedirs(args.out_dir_p2, exist_ok=True)

    input_files = sorted(glob.glob(os.path.join(args.input_dir, "*.txt")))

    if not input_files:
        print(f"Aucun fichier .txt trouve dans {args.input_dir}")
        return

    for input_path in input_files:
        filename = os.path.basename(input_path)

        out_3 = os.path.join(args.out_dir_p3, filename)
        out_2 = os.path.join(args.out_dir_p2, filename)

        truncate_instance(input_path, out_3, remove_last_n=1)
        truncate_instance(input_path, out_2, remove_last_n=2)

        print(f"{filename} -> {out_3} (p-1 criteres) et {out_2} (p-2 criteres)")

    print(f"\nTermine : {len(input_files)} instance(s) traitee(s).")


if __name__ == "__main__":
    main()
