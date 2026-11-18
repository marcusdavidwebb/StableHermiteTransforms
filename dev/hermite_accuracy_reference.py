"""Generate high-precision reference factors for the Hermite accuracy comparison.

Usage: python3 hermite_accuracy_reference.py input_nodes.txt output.bin
Input is one approximate Gauss-Hermite node per line, in increasing order.
Output contains the nonnegative nodes' records, in increasing node order.
Each record is N+1 little-endian binary64 numbers: d followed by Q[:, node].
All calculations before writing use 64 decimal digits.
"""

import argparse
import sys
from array import array
from decimal import Decimal, localcontext
from pathlib import Path


PI = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406")


def refine_root(seed, half_degrees):
    """Return a root and its monic Hermite polynomial values up to degree N."""
    n = len(half_degrees) + 1
    x = seed
    for _ in range(8):
        # p_(k+1) = x p_k - (k/2) p_(k-1), and p'_N = N p_(N-1).
        polynomials = [Decimal(1), x]
        previous, current = polynomials
        for half_degree in half_degrees:
            previous, current = current, x * current - half_degree * previous
            polynomials.append(current)
        correction = current / (n * previous)
        if abs(correction) <= Decimal("1e-54") * max(1, abs(x)):
            return x, polynomials
        x -= correction
    raise ArithmeticError(f"Root refinement failed for N={n}, seed={seed}")


def write_reference(input_path, output_path):
    seeds = [Decimal(line) for line in input_path.read_text().splitlines()]
    n = len(seeds)
    if not seeds or not all(x.is_finite() for x in seeds):
        raise ValueError("Expected a nonempty list of finite Gauss-Hermite nodes")
    if any(left >= right for left, right in zip(seeds, seeds[1:])):
        raise ValueError("Expected Gauss-Hermite nodes in strictly increasing order")

    with localcontext() as context, output_path.open("wb") as output:
        context.prec = 64
        half_degrees = [Decimal(k) / 2 for k in range(1, n)]
        sqrt_n = Decimal(n).sqrt()
        # The Hermite function of degree k is exp(-x^2/2) * c_k * p_k(x).
        coefficients = [1 / PI.sqrt().sqrt()]
        for k in range(1, n):
            coefficients.append(coefficients[-1] * (Decimal(2) / k).sqrt())

        positive_half = seeds[n // 2 :]
        if n % 2:
            positive_half[0] = Decimal(0)
        last_root = None
        for seed in positive_half:
            root, polynomials = refine_root(seed, half_degrees)
            if root < 0 or (last_root is not None and root <= last_root):
                raise ArithmeticError(f"Expected distinct nonnegative roots at N={n}")
            last_root = root

            # The common Gaussian cancels in Q. Decimal's exponent range
            # accommodates the unscaled polynomials at every size in the paper.
            divisor = sqrt_n * abs(coefficients[-1] * polynomials[n - 1])
            d = (-root * root / 2).exp() * divisor
            record = array("d", [float(d)])
            record.extend(float(c * p / divisor) for c, p in zip(coefficients, polynomials))
            if sys.byteorder != "little":
                record.byteswap()
            record.tofile(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_nodes", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    write_reference(args.input_nodes, args.output)
