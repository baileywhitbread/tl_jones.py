"""Finite-field reconstruction of Temperley-Lieb Jones traces.

This module complements tl_jones.py by reconstructing exact Laurent
polynomials from finite-field evaluations of dense braid matrices.
"""

from __future__ import annotations
from dataclasses import dataclass
from functools import cache
from typing import Iterable, Sequence

import sympy as sp
from sympy.ntheory import isprime, primitive_root

from tl_jones import (
    _apply_temperley_lieb_generator,
    _validate_word,
    link_states,
    quantum_integer,
    quantum_two,
    v,
)

FFMatrix = list[list[int]]

@dataclass(frozen=True)
class FFEvaluationContext:
    """Shared context for reconstructing the trace of a braid in a Jones representation."""
    n: int
    r: int
    word: tuple[int, ...]
    p: int
    order: int
    root_of_unity: int
    degree_bound: int
    coefficient_bound: int
    dimension: int


# Finite-field analogues of functions in tl_jones.py.
def ff_jones_polynomial(n: int, word: Iterable[int]) -> sp.Expr:
    """Return the unnormalised Jones polynomial via finite-field traces."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    word = _validate_word(n, word)
    writhe = sum(1 if generator > 0 else -1 for generator in word)
    total = sp.Integer(0)
    for r in range(n // 2 + 1):
        coefficient = quantum_integer(n - 2 * r + 1)
        total += coefficient * ff_jones_rep_trace(n, r, word)

    return sp.expand(((-1) ** writhe) * v ** (-2 * writhe) * total)


def ff_jones_rep_trace(n: int, r: int, word: Iterable[int]) -> sp.Expr:
    """Return tr(rho_(n-r,r)(word)) by finite-field interpolation."""
    context = _ff_create_context(n, r, word)
    braid_mats = ff_jones_rep_braid_word(n, r, context.word)
    traces = [
        sum(matrix[index][index] for index in range(context.dimension)) % context.p
        for matrix in braid_mats
    ]

    residues = _inverse_dft(
        traces,
        context.p,
        context.root_of_unity,
        context.degree_bound,
    )
    expression = sp.Integer(0)
    for degree, residue in enumerate(residues, start=-context.degree_bound):
        residue %= context.p
        if residue > context.p // 2:
            lifted = residue - context.p
        else:
            lifted = residue
        if lifted != 0:
            expression += sp.Integer(lifted) * v**degree
    return sp.expand(expression)

def ff_jones_rep_braid_word(n: int, r: int, word: Iterable[int]) -> list[FFMatrix]:
    """Evaluate and multiply a braid word at every interpolation root."""
    context = _ff_create_context(n, r, word)
    action_by_index = {
        index: _ff_temperley_lieb_matrix(context.n, context.r, index)
        for index in {abs(generator) for generator in context.word}
    }

    braid_mats: list[FFMatrix] = []
    current_root = 1
    for _ in range(context.order):
        inverse_root = pow(current_root, -1, context.p)
        delta = (current_root + inverse_root) % context.p

        matrix = _identity_matrix(context.dimension)
        for generator in context.word:
            scalar = inverse_root if generator < 0 else current_root
            generator_matrix = [
                [0] * context.dimension for _ in range(context.dimension)
            ]

            for column, (kind, row) in enumerate(action_by_index[abs(generator)]):
                generator_matrix[column][column] = (-scalar) % context.p
                if kind == 1:
                    generator_matrix[row][column] = (
                        generator_matrix[row][column] + 1
                    ) % context.p
                elif kind == quantum_two:
                    generator_matrix[row][column] = (
                        generator_matrix[row][column] + delta
                    ) % context.p

            matrix = _matrix_multiply_mod(matrix, generator_matrix, context.p)
        braid_mats.append(matrix)
        current_root = (current_root * context.root_of_unity) % context.p

    return braid_mats

# Auxiliary helpers.
def _ff_create_context(n: int, r: int, word: Iterable[int]) -> FFEvaluationContext:
    """Build the finite-field reconstruction context for one trace."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    if type(r) is not int or r < 0 or 2 * r > n:
        raise ValueError("r must be an integer with 0 <= r <= floor(n / 2)")
    word = _validate_word(n, word)

    dimension = len(link_states(n, r))

    if len(word) == 0:
        degree_bound = 0
        order = 1
        coefficient_bound = dimension
        p = int(sp.nextprime(2 * coefficient_bound))
        root_of_unity = 1
    else:
        degree_bound = len(word)
        order = 2 * degree_bound + 1
        coefficient_bound = dimension ** len(word)
        candidate = ((2 * coefficient_bound) // order) * order + 1
        if candidate <= (2 * coefficient_bound):
            candidate += order
        while not isprime(candidate):
            candidate += order
        p = candidate

        generator = primitive_root(p)
        root_of_unity = pow(int(generator), (p - 1) // order, p)

    return FFEvaluationContext(
        n=n,
        r=r,
        word=word,
        p=p,
        order=order,
        root_of_unity=root_of_unity,
        degree_bound=degree_bound,
        coefficient_bound=coefficient_bound,
        dimension=dimension,
    )


def _inverse_dft(
    values: Sequence[int], p: int, root_of_unity: int, bound: int
) -> list[int]:
    """Recover Laurent coefficients modulo p from root evaluations."""
    order = len(values)
    shifted_values = [
        (value * pow(root_of_unity, index * bound, p)) % p
        for index, value in enumerate(values)
    ]

    inverse_root = pow(root_of_unity, -1, p)
    inverse_order = pow(order, -1, p)
    coefficients = [0] * order

    for degree in range(order):
        base = pow(inverse_root, degree, p)
        power = 1
        total = 0
        for value in shifted_values:
            total = (total + value * power) % p
            power = (power * base) % p
        coefficients[degree] = (total * inverse_order) % p

    return coefficients

@cache
def _ff_temperley_lieb_matrix(
    n: int, r: int, i: int
) -> tuple[tuple[sp.Expr, int], ...]:
    """
    Return the sparse column action of E_i on the link-state basis.
    Note: We are reimplementing the same logic as tl_jones._temperley_lieb_matrix
    but SymPy expressions are avoided and the result is cached for efficiency.
    """

    if type(i) is not int or i < 1 or i >= n:
        raise ValueError("generator index i must satisfy 1 <= i < n")

    basis = link_states(n, r)
    state_to_index = {state: index for index, state in enumerate(basis)}
    action_data: list[tuple[sp.Expr, int]] = []

    for state in basis:
        coefficient, next_state = _apply_temperley_lieb_generator(state, i)
        if coefficient == 0:
            action_data.append((coefficient, -1))
        else:
            action_data.append((coefficient, state_to_index[next_state]))

    return tuple(action_data)



def _identity_matrix(size: int) -> FFMatrix:
    """Return the size-by-size identity matrix as nested lists."""
    matrix = [[0] * size for _ in range(size)]
    for index in range(size):
        matrix[index][index] = 1
    return matrix


def _matrix_multiply_mod(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]], p: int
) -> FFMatrix:
    """Multiply two dense matrices modulo p."""
    size = len(left)
    result = [[0] * size for _ in range(size)]

    for row_index, left_row in enumerate(left):
        result_row = result[row_index]
        for middle_index, left_entry in enumerate(left_row):
            if left_entry == 0:
                continue
            right_row = right[middle_index]
            for column_index, right_entry in enumerate(right_row):
                if right_entry == 0:
                    continue
                result_row[column_index] = (
                    result_row[column_index] + left_entry * right_entry
                ) % p

    return result


__all__ = [
    "FFMatrix",
    "FFEvaluationContext",
    "ff_jones_polynomial",
    "ff_jones_rep_trace",
    "ff_jones_rep_braid_word",
]
