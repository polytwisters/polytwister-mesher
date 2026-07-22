import numpy as np
import json
import argparse


def similarity(a, b):
    return np.abs(np.dot(a, b.conj())) / (np.linalg.norm(a) * np.linalg.norm(b))


def random_c2(rng):
    result = rng.standard_normal(4)
    result = result / np.linalg.norm(result)
    radius = rng.triangular(0.7, 1.0, 1.3)
    result = result * radius
    result = np.round(result, 1)
    return r4_to_c2(result)

def r4_to_c2(r4):
    return r4[0::2] + 1j * r4[1::2]

def c2_to_r4(c2):
    result = np.zeros(len(c2) * 2)
    result[0::2] = c2.real
    result[1::2] = c2.imag
    return result


def random_convex_polytwister(n: int):
    logs_c2 = []
    logs_real = []
    rng = np.random.default_rng(0)
    for i in range(n):
        c2 = random_c2(rng)
        for __ in range(100):
            max_similarity = max([similarity(c2, log_c2) for log_c2 in logs_c2], default=0.0)
            if max_similarity < 0.9:
                break
        logs_c2.append(c2)
        logs_real.append(c2_to_r4(c2).tolist())
    return {"logs": logs_real}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int)
    args = parser.parse_args()

    print(json.dumps(random_convex_polytwister(args.n)))


if __name__ == "__main__":
    main()