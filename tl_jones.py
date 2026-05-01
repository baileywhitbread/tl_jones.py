"""Jones polynomials from braid words via Temperley-Lieb representations.

This module follows the setup in docs/documentation.tex.
It computes the unnormalised Jones polynomial of a braid closure using the
link-state representations of the Temperley-Lieb algebra.
"""

from __future__ import annotations
from functools import cache
from typing import Iterable
import sympy as sp

v = sp.Symbol("v")
quantum_two = v + v**-1


def quantum_integer(m: int) -> sp.Expr:
    """Return the quantum integer [m]_v."""
    if type(m) is not int or m < 0:
        raise ValueError("m must be a nonnegative integer")
    return sp.expand(sum(v ** (m - 1 - 2 * k) for k in range(m)))


def projlen(polynomial: sp.Expr) -> int:
    """Return the spread between the minimum and maximum v-degrees."""
    polynomial = sp.expand(polynomial)
    if polynomial == 0:
        return 0

    degrees = [
        int(term.as_powers_dict().get(v, 0))
        for term in polynomial.as_ordered_terms()
    ]
    return max(degrees) - min(degrees)


def gnf_to_braid_word(
    n: int, d: int, factors: Iterable[Iterable[int]]
) -> list[int]:
    """Expand classical Garside normal form data to a braid word."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    if type(d) is not int:
        raise ValueError("d must be an integer")
    factors = _validate_gnf(n, factors)

    garside_word = _garside_delta_word(n)
    if d < 0:
        garside_word = [-generator for generator in reversed(garside_word)]

    braid_word = garside_word * abs(d)
    for factor in factors:
        braid_word += _permutation_to_braid_word(factor)

    return braid_word


@cache
def link_states(n: int, r: int) -> tuple[tuple[int, ...], ...]:
    """Return the link-state basis for V_(n-r,r).

    A state is a tuple of length n.  The entry -1 marks a defect, and a
    nonnegative entry records the index of the paired endpoint.
    """
    if type(n) is not int or n < 0:
        raise ValueError("n must be a nonnegative integer")
    if type(r) is not int or r < 0 or 2 * r > n:
        raise ValueError("r must be an integer with 0 <= r <= floor(n / 2)")

    if n == 0:
        return ((),)

    states: list[tuple[int, ...]] = []

    # First allow the leftmost point to be a defect, so the basis is ordered
    # with defects as far left as possible.
    if 2 * r <= n - 1:
        for rest in link_states(n - 1, r):
            shifted_rest = tuple(
                -1 if entry == -1 else entry + 1 for entry in rest
            )
            states.append((-1,) + shifted_rest)

    # Then make the leftmost point part of a cap. The points inside that cap
    # must be fully paired, otherwise a defect would cross the cap.
    for partner in range(1, n):
        inside_length = partner - 1
        if inside_length % 2 != 0:
            continue

        inside_caps = inside_length // 2
        outside_length = n - partner - 1
        outside_caps = r - inside_caps - 1
        if outside_caps < 0 or 2 * outside_caps > outside_length:
            continue

        for inside in link_states(inside_length, inside_caps):
            for outside in link_states(outside_length, outside_caps):
                state = [None] * n
                state[0] = partner
                state[partner] = 0

                for local_index, paired_with in enumerate(inside, start=1):
                    state[local_index] = (
                        -1 if paired_with == -1 else paired_with + 1
                    )

                outside_offset = partner + 1
                for local_index, paired_with in enumerate(
                    outside, start=outside_offset
                ):
                    state[local_index] = (
                        -1 if paired_with == -1
                        else paired_with + outside_offset
                    )

                states.append(tuple(state))

    return tuple(states)


def jones_rep_braid_generator(
    n: int, r: int, i: int, inverse: bool = False
) -> sp.MatrixBase:
    """Return rho_(n-r,r)(sigma_i^(+/-1)) via E_i - v^(+/-1) I."""
    e_i = _temperley_lieb_matrix(n, r, i)
    identity = sp.eye(e_i.rows)
    scalar = v**-1 if inverse else v
    return e_i - scalar * identity


def jones_rep_braid_word(n: int, r: int, word: Iterable[int]) -> sp.Matrix:
    """Return rho_(n-r,r)(word) for a signed-integer braid word.

    Positive entries are positive braid generators.  Negative entries are
    inverse braid generators.  For example, [1, -2, 1] means
    sigma_1 sigma_2^-1 sigma_1.
    """
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    if type(r) is not int or r < 0 or 2 * r > n:
        raise ValueError("r must be an integer with 0 <= r <= floor(n / 2)")
    word = _validate_word(n, word)

    dimension = len(link_states(n, r))
    matrix = sp.eye(dimension)

    for generator in word:
        index = abs(generator)
        factor = jones_rep_braid_generator(
            n, r, index, inverse=(generator < 0)
        )
        matrix = matrix * factor

    return matrix


def jones_polynomial(n: int, word: Iterable[int]) -> sp.Expr:
    """Return the unnormalised Jones polynomial via a writhe-corrected trace."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    word = _validate_word(n, word)
    writhe = sum(1 if generator > 0 else -1 for generator in word)

    total = sp.Integer(0)
    for r in range(n // 2 + 1):
        coefficient = quantum_integer(n - 2 * r + 1)
        total += coefficient * sp.trace(jones_rep_braid_word(n, r, word))

    return sp.expand(((-1) ** writhe) * v ** (-2 * writhe) * total)


@cache
def _temperley_lieb_matrix(n: int, r: int, i: int) -> sp.ImmutableMatrix:
    """Return the matrix of E_i on the link-state basis for V_(n-r,r)."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    if type(r) is not int or r < 0 or 2 * r > n:
        raise ValueError("r must be an integer with 0 <= r <= floor(n / 2)")
    if type(i) is not int or i < 1 or i >= n:
        raise ValueError("generator index i must satisfy 1 <= i < n")

    basis = link_states(n, r)
    state_to_index = {state: index for index, state in enumerate(basis)}
    matrix = sp.zeros(len(basis))

    for column, state in enumerate(basis):
        coefficient, next_state = _apply_temperley_lieb_generator(state, i)
        if coefficient != 0:
            row = state_to_index[next_state]
            matrix[row, column] += coefficient

    # Freeze the matrix before caching so later callers cannot mutate it.
    return sp.ImmutableMatrix(matrix)


def _apply_temperley_lieb_generator(
    state: tuple[int, ...], i: int
) -> tuple[sp.Expr, tuple[int, ...]]:
    """Apply the local Temperley-Lieb generator E_i to one link state."""
    left = i - 1
    right = i
    left_partner = state[left]
    right_partner = state[right]

    # If there is a cap at (i, i+1), we equal the same state with a loop.
    if left_partner == right and right_partner == left:
        return quantum_two, state

    # If there are defects at i and i+1, we have left the link-state basis.
    if left_partner == -1 and right_partner == -1:
        return sp.Integer(0), state

    # Otherwise, no loops are created, and we have not left the link-state
    # basis, so we must equal a single link-state basis element.
    new_state = list(state)
    new_state[left] = right
    new_state[right] = left

    # If exactly one point was a defect, the far endpoint of the old pair
    # becomes the new defect.
    if left_partner == -1:
        new_state[right_partner] = -1
    elif right_partner == -1:
        new_state[left_partner] = -1
    else:
        # If both points were paired elsewhere, reconnect the two displaced
        # endpoints to each other.
        new_state[left_partner] = right_partner
        new_state[right_partner] = left_partner

    return sp.Integer(1), tuple(new_state)


def _garside_delta_word(n: int) -> list[int]:
    """Return the chosen braid word for the Garside element Delta_n."""
    word = []
    for last_generator in range(n - 1, 0, -1):
        word += list(range(1, last_generator + 1))
    return word


def _permutation_to_braid_word(permutation: tuple[int, ...]) -> list[int]:
    """Convert a one-line permutation into a positive braid word."""
    current = list(range(1, len(permutation) + 1))
    word = []

    for target_index, desired_value in enumerate(permutation):
        current_index = current.index(desired_value)
        while current_index > target_index:
            word.append(current_index)
            current[current_index - 1], current[current_index] = (
                current[current_index],
                current[current_index - 1],
            )
            current_index -= 1

    return word


def _validate_gnf(
    n: int, factors: Iterable[Iterable[int]]
) -> tuple[tuple[int, ...], ...]:
    """Check the full Garside normal-form input for gnf_to_braid_word.

    This verifies that each factor is a one-line permutation of 1, ..., n and
    rules out the identity and Delta permutations. It also checks that
    consecutive simple factors satisfy the Garside descent-set rule before the
    GNF is expanded to a braid word.
    """
    try:
        factors = tuple(tuple(factor) for factor in factors)
    except TypeError as error:
        raise ValueError(
            "factors must be an iterable of one-line permutations"
        ) from error

    identity = tuple(range(1, n + 1))
    delta_permutation = tuple(range(n, 0, -1))

    for factor in factors:
        if len(factor) != n:
            raise ValueError("each factor must be a permutation of 1, ..., n")
        if any(type(entry) is not int for entry in factor):
            raise ValueError("each factor must be a permutation of 1, ..., n")
        if sorted(factor) != list(range(1, n + 1)):
            raise ValueError("each factor must be a permutation of 1, ..., n")
        if factor == identity:
            raise ValueError("identity is not allowed as a simple factor")
        if factor == delta_permutation:
            raise ValueError("Delta is not allowed as a simple factor")

    for previous, next_factor in zip(factors, factors[1:]):
        if not _left_descent_set(next_factor).issubset(
            _right_descent_set(previous)
        ):
            raise ValueError("factors do not satisfy the Garside descent-set rule")

    return factors


def _inverse_permutation(permutation: tuple[int, ...]) -> tuple[int, ...]:
    """Return the inverse of a one-line permutation."""
    inverse = [0] * len(permutation)
    for index, value in enumerate(permutation, start=1):
        inverse[value - 1] = index
    return tuple(inverse)


def _right_descent_set(permutation: tuple[int, ...]) -> set[int]:
    """Return the right descent set of a one-line permutation."""
    return {
        index
        for index in range(1, len(permutation))
        if permutation[index - 1] > permutation[index]
    }


def _left_descent_set(permutation: tuple[int, ...]) -> set[int]:
    """Return the left descent set of a one-line permutation."""
    return _right_descent_set(_inverse_permutation(permutation))


def _validate_word(n: int, word: Iterable[int]) -> tuple[int, ...]:
    """Check that word is a braid word in nonzero integer generators.

    This is used by jones_rep_braid_word and jones_polynomial, and it
    also checks that each generator index lies in the range 1, ..., n - 1.
    """
    try:
        word = tuple(word)
    except TypeError as error:
        raise ValueError("word must be an iterable of nonzero integers") from error

    for generator in word:
        if type(generator) is not int:
            raise ValueError("braid word entries must be nonzero integers")
        if generator == 0:
            raise ValueError("braid word entries must be nonzero integers")
        if abs(generator) >= n:
            raise ValueError("generator index i must satisfy 1 <= i < n")

    return word


__all__ = [
    "v",
    "quantum_two",
    "quantum_integer",
    "projlen",
    "gnf_to_braid_word",
    "link_states",
    "jones_rep_braid_generator",
    "jones_rep_braid_word",
    "jones_polynomial",
]
