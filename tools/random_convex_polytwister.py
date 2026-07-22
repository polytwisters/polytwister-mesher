import numpy as np
import json
import argparse


def similarity(a, b):
    return np.abs(np.dot(a, b.conj())) / (np.linalg.norm(a) * np.linalg.norm(b))


def random_vec4(rng):
    r4 = rng.standard_normal(4)
    r4 = r4 / np.linalg.norm(r4)
    radius = rng.triangular(0.8, 1.0, 1.2)
    r4 = r4 * radius
    r4 = np.round(r4, 1)
    c2 = r4[0::2] + 1j * r4[1::2]
    return (c2, r4)


def random_convex_polytwister(n: int):
    logs_c2 = []
    logs_real = []
    rng = np.random.default_rng(0)
    for i in range(n):
        c2, r4 = random_vec4(rng)
        for __ in range(1000):
            if min([similarity(c2, log_c2) for log_c2 in logs_c2], default=0.0) < 0.9:
                break
            c2, r4 = random_vec4(rng)
        else:
            print("warning not found")
        logs_c2.append(c2)
        logs_real.append(r4.tolist())
    return {"logs": logs_real}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int)
    args = parser.parse_args()

    print(json.dumps(random_convex_polytwister(args.n)))


if __name__ == "__main__":
    main()