import ast
import csv
import random
import sys
import unittest
from pathlib import Path

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import tl_jones as tl
import tl_jones_ff as tlff

FIXTURE_PATH = PROJECT_ROOT / "data" / "jones_braid_knots.csv"


class FiniteFieldTraceTests(unittest.TestCase):
    def assertExprEqual(self, left, right):
        self.assertEqual(sp.simplify(left - right), 0)

    def assertMatrixShape(self, matrix, rows, columns):
        self.assertEqual(len(matrix), rows)
        for row in matrix:
            self.assertEqual(len(row), columns)

    def assertMatricesEqual(self, left, right):
        self.assertEqual(len(left), len(right))
        for left_matrix, right_matrix in zip(left, right):
            self.assertEqual(left_matrix, right_matrix)

    def random_braid_word(self, rng, n, length):
        if n == 1:
            return []
        return [
            rng.choice([-1, 1]) * rng.randrange(1, n)
            for _ in range(length)
        ]

    def load_knotinfo_fixture(self):
        with FIXTURE_PATH.open(newline="") as file:
            rows = []
            for row in csv.DictReader(file):
                rows.append(
                    {
                        "name": row["name"],
                        "crossing_number": int(row["crossing_number"]),
                        "strands": int(row["strands"]),
                        "braid_word": ast.literal_eval(row["braid_word"]),
                        "jones_q": row["jones_q"],
                    }
                )
            return rows

    def test_01_validation_errors(self):
        with self.assertRaises(ValueError):
            tlff.ff_jones_rep_trace(0, 0, [])
        with self.assertRaises(ValueError):
            tlff.ff_jones_rep_trace(3, 2, [])
        with self.assertRaises(ValueError):
            tlff.ff_jones_rep_trace(2, 1, [0])
        with self.assertRaises(ValueError):
            tlff.ff_jones_rep_trace(2, 1, [2])
        with self.assertRaises(ValueError):
            tlff.ff_jones_rep_trace(2, 1, [1.0])
        with self.assertRaises(ValueError):
            tlff.ff_jones_polynomial(0, [])
        with self.assertRaises(ValueError):
            tlff.ff_jones_polynomial(2, [0])
        with self.assertRaises(ValueError):
            tlff.ff_jones_rep_braid_word(2, 1, [1.0])

    def test_02_reconstruction_prime_is_large_enough_and_has_required_order(self):
        cases = [
            (5, 2, []),
            (2, 1, [1, 1, 1]),
            (3, 1, [1, -2, 1, -2]),
            (4, 2, [1, 2, -1, 3, -2, 1]),
        ]

        for n, r, word in cases:
            with self.subTest(n=n, r=r, word=word):
                dimension = len(tl.link_states(n, r))
                degree_bound = len(word)
                order = 1 if len(word) == 0 else 2 * degree_bound + 1
                coefficient_bound = (
                    dimension if len(word) == 0 else dimension ** len(word)
                )
                context = tlff._ff_create_context(n, r, word)
                self.assertEqual(context.n, n)
                self.assertEqual(context.r, r)
                self.assertEqual(context.word, tuple(word))
                self.assertEqual(context.dimension, dimension)
                self.assertEqual(context.degree_bound, degree_bound)
                self.assertEqual(context.order, order)
                self.assertEqual(context.coefficient_bound, coefficient_bound)
                self.assertTrue(sp.isprime(context.p))
                self.assertGreater(context.p, 2 * coefficient_bound)
                self.assertGreaterEqual(context.root_of_unity, 0)
                self.assertLess(context.root_of_unity, context.p)
                self.assertEqual(
                    pow(context.root_of_unity, context.order, context.p),
                    1,
                )
                if context.order > 1:
                    self.assertEqual(context.p % context.order, 1)
                    for exponent in range(1, context.order):
                        self.assertNotEqual(
                            pow(context.root_of_unity, exponent, context.p),
                            1,
                        )

    def test_02a_braid_matrix_stage_shapes(self):
        n, r, word = 4, 1, [1, 2, -1, 3]
        context = tlff._ff_create_context(n, r, word)
        braid_mats = tlff.ff_jones_rep_braid_word(n, r, word)

        self.assertEqual(context.n, n)
        self.assertEqual(context.r, r)
        self.assertEqual(context.word, tuple(word))
        self.assertEqual(context.degree_bound, len(word))
        self.assertEqual(context.order, 2 * len(word) + 1)
        self.assertEqual(len(braid_mats), context.order)

        for matrix in braid_mats:
            self.assertMatrixShape(
                matrix, context.dimension, context.dimension
            )
            for row in matrix:
                for entry in row:
                    self.assertIs(type(entry), int)
                    self.assertGreaterEqual(entry, 0)
                    self.assertLess(entry, context.p)

    def test_02b_four_stage_pipeline_matches_trace_wrapper(self):
        cases = [
            (2, 1, [1, 1, 1]),
            (3, 1, [1, -2, 1, -2]),
            (4, 2, [1, 2, -1, 3, -2, 1]),
        ]

        for n, r, word in cases:
            with self.subTest(n=n, r=r, word=word):
                context = tlff._ff_create_context(n, r, word)
                braid_mats = tlff.ff_jones_rep_braid_word(n, r, word)
                traces = [
                    sum(matrix[index][index] for index in range(context.dimension))
                    % context.p
                    for matrix in braid_mats
                ]
                residues = tlff._inverse_dft(
                    traces,
                    context.p,
                    context.root_of_unity,
                    context.degree_bound,
                )
                actual = sp.Integer(0)
                for degree, residue in enumerate(
                    residues, start=-context.degree_bound
                ):
                    residue %= context.p
                    lifted = (
                        residue - context.p
                        if residue > context.p // 2
                        else residue
                    )
                    if lifted != 0:
                        actual += sp.Integer(lifted) * tl.v**degree
                actual = sp.expand(actual)
                expected = tlff.ff_jones_rep_trace(n, r, word)
                self.assertExprEqual(actual, expected)

    def test_02c_empty_word_flows_through_all_pipeline_stages(self):
        n, r, word = 5, 2, []
        dimension = len(tl.link_states(n, r))

        context = tlff._ff_create_context(n, r, word)
        braid_mats = tlff.ff_jones_rep_braid_word(n, r, word)
        self.assertEqual(context.degree_bound, 0)
        self.assertEqual(context.order, 1)
        self.assertEqual(context.root_of_unity, 1)
        self.assertEqual(context.coefficient_bound, dimension)
        self.assertEqual(braid_mats, [tlff._identity_matrix(dimension)])

        traces = [
            sum(matrix[index][index] for index in range(context.dimension))
            % context.p
            for matrix in braid_mats
        ]
        self.assertEqual(traces, [dimension])

        residues = tlff._inverse_dft(
            traces,
            context.p,
            context.root_of_unity,
            context.degree_bound,
        )
        actual = sp.Integer(0)
        for degree, residue in enumerate(residues, start=-context.degree_bound):
            residue %= context.p
            lifted = residue - context.p if residue > context.p // 2 else residue
            if lifted != 0:
                actual += sp.Integer(lifted) * tl.v**degree
        actual = sp.expand(actual)
        self.assertEqual(actual, sp.Integer(dimension))

    def test_03_ff_jones_polynomial_matches_symbolic_backend_small_cases(self):
        cases = [
            (2, [1, 1, 1]),
            (3, [1, -2, 1, -2]),
            (4, [1, 2, -1, 3, -2, 1]),
        ]

        for n, word in cases:
            with self.subTest(n=n, word=word):
                expected = tl.jones_polynomial(n, word)
                actual = tlff.ff_jones_polynomial(n, word)
                self.assertExprEqual(actual, expected)

    def test_03a_ff_jones_polynomial_matches_symbolic_backend_random_cases(self):
        rng = random.Random(20260430)

        for case_index in range(12):
            n = rng.randrange(1, 6)
            word = self.random_braid_word(rng, n, rng.randrange(0, 6))
            with self.subTest(case_index=case_index, n=n, word=word):
                expected = tl.jones_polynomial(n, word)
                actual = tlff.ff_jones_polynomial(n, word)
                self.assertExprEqual(actual, expected)

    def test_04_inverse_dft_recovers_random_laurent_polynomials(self):
        rng = random.Random(20260430)

        for bound in range(5):
            order = 2 * bound + 1
            if order == 1:
                p = 101
                root_of_unity = 1
            else:
                p = 101
                while p % order != 1:
                    p = int(sp.nextprime(p))
                root_of_unity = pow(
                    int(sp.primitive_root(p)),
                    (p - 1) // order,
                    p,
                )

            for trial in range(5):
                coefficients = [
                    rng.randrange(-25, 26)
                    for _ in range(order)
                ]
                values = []
                for sample_index in range(order):
                    point = pow(root_of_unity, sample_index, p)
                    value = 0
                    for degree, coefficient in enumerate(
                        coefficients, start=-bound
                    ):
                        value += coefficient * pow(point, degree, p)
                    values.append(value % p)

                with self.subTest(bound=bound, trial=trial):
                    residues = tlff._inverse_dft(
                        values,
                        p,
                        root_of_unity,
                        bound,
                    )
                    self.assertEqual(
                        residues,
                        [coefficient % p for coefficient in coefficients],
                    )

    def test_05_ff_jones_rep_trace_matches_symbolic_backend_small_cases(self):
        cases = [
            (2, 0, [1, 1, 1]),
            (2, 1, [1, 1, 1]),
            (3, 1, [1, -2, 1, -2]),
            (3, 1, [1, 1, 1, 2, -1, 2]),
            (4, 1, [1, 2, -1, 3, -2, 1]),
        ]

        for n, r, word in cases:
            with self.subTest(n=n, r=r, word=word):
                expected = sp.expand(sp.trace(tl.jones_rep_braid_word(n, r, word)))
                actual = tlff.ff_jones_rep_trace(n, r, word)
                self.assertExprEqual(actual, expected)

    def test_05b_ff_jones_rep_trace_matches_symbolic_backend_random_cases(self):
        rng = random.Random(20260430)

        for n in range(1, 6):
            for case_index in range(3):
                word = self.random_braid_word(rng, n, rng.randrange(0, 6))
                for r in range(n // 2 + 1):
                    with self.subTest(
                        n=n,
                        r=r,
                        case_index=case_index,
                        word=word,
                    ):
                        expected = sp.expand(
                            sp.trace(tl.jones_rep_braid_word(n, r, word))
                        )
                        actual = tlff.ff_jones_rep_trace(n, r, word)
                        self.assertExprEqual(actual, expected)

    def test_05a_ff_jones_rep_trace_can_have_positive_and_negative_degrees(self):
        actual = sp.expand(tlff.ff_jones_rep_trace(3, 1, [1, -2, 1, -2]))
        degrees = [
            int(term.as_powers_dict().get(tl.v, 0))
            for term in actual.as_ordered_terms()
        ]
        self.assertLess(min(degrees), 0)
        self.assertGreater(max(degrees), 0)

    def test_06_braid_inverse_relation(self):
        for n in range(2, 8):
            for r in range(n // 2 + 1):
                dimension = sp.Integer(len(tl.link_states(n, r)))
                for i in range(1, n):
                    with self.subTest(n=n, r=r, i=i, word=[i, -i]):
                        self.assertExprEqual(
                            tlff.ff_jones_rep_trace(n, r, [i, -i]),
                            dimension,
                        )

                    with self.subTest(n=n, r=r, i=i, word=[-i, i]):
                        self.assertExprEqual(
                            tlff.ff_jones_rep_trace(n, r, [-i, i]),
                            dimension,
                        )

    def test_07_braid_commuting_relation(self):
        for n in range(4, 8):
            for r in range(n // 2 + 1):
                for i in range(1, n):
                    for j in range(1, n):
                        if abs(i - j) > 1:
                            left = tlff.ff_jones_rep_trace(n, r, [i, j])
                            right = tlff.ff_jones_rep_trace(n, r, [j, i])
                            with self.subTest(n=n, r=r, i=i, j=j):
                                self.assertExprEqual(left, right)

    def test_08_braid_relation(self):
        for n in range(3, 8):
            for r in range(n // 2 + 1):
                for i in range(1, n - 1):
                    left = tlff.ff_jones_rep_trace(n, r, [i, i + 1, i])
                    right = tlff.ff_jones_rep_trace(
                        n, r, [i + 1, i, i + 1]
                    )
                    with self.subTest(n=n, r=r, i=i):
                        self.assertExprEqual(left, right)

    def test_08a_evaluated_braid_matrices_satisfy_braid_relations(self):
        for n in range(3, 6):
            for r in range(n // 2 + 1):
                for i in range(1, n):
                    context = tlff._ff_create_context(n, r, [i, -i])
                    identity = tlff._identity_matrix(context.dimension)
                    matrices = tlff.ff_jones_rep_braid_word(n, r, [i, -i])
                    with self.subTest(relation="inverse", n=n, r=r, i=i):
                        self.assertEqual(
                            matrices,
                            [identity for _ in range(context.order)],
                        )

                for i in range(1, n - 1):
                    with self.subTest(relation="yang-baxter", n=n, r=r, i=i):
                        self.assertMatricesEqual(
                            tlff.ff_jones_rep_braid_word(n, r, [i, i + 1, i]),
                            tlff.ff_jones_rep_braid_word(n, r, [i + 1, i, i + 1]),
                        )

                for i in range(1, n):
                    for j in range(1, n):
                        if abs(i - j) > 1:
                            with self.subTest(
                                relation="commuting",
                                n=n,
                                r=r,
                                i=i,
                                j=j,
                            ):
                                self.assertMatricesEqual(
                                    tlff.ff_jones_rep_braid_word(n, r, [i, j]),
                                    tlff.ff_jones_rep_braid_word(n, r, [j, i]),
                                )

    def test_09_ff_jones_rep_trace_simple_cases(self):
        self.assertEqual(tlff.ff_jones_rep_trace(1, 0, []), 1)

        for n in range(1, 11):
            for r in range(n // 2 + 1):
                with self.subTest(n=n, r=r, word=[]):
                    self.assertEqual(
                        tlff.ff_jones_rep_trace(n, r, []),
                        sp.Integer(len(tl.link_states(n, r))),
                    )

        cases = [
            ((2, 0, [1]), -tl.v),
            ((2, 1, [1]), tl.v**-1),
            ((2, 0, [-1]), -tl.v**-1),
            ((2, 1, [-1]), tl.v),
            ((2, 0, [1, 1]), tl.v**2),
            ((2, 1, [1, 1]), tl.v**-2),
            ((2, 0, [1, 1, 1]), -tl.v**3),
            ((2, 1, [1, 1, 1]), tl.v**-3),
            ((3, 1, [1]), -tl.v + tl.v**-1),
            ((3, 1, [-1]), tl.v - tl.v**-1),
        ]

        for (n, r, word), expected in cases:
            with self.subTest(n=n, r=r, word=word):
                actual = tlff.ff_jones_rep_trace(n, r, word)
                self.assertExprEqual(actual, expected)

    def test_09a_one_shot_iterables_are_accepted(self):
        word = [1, -2, 1, 2]

        trace_from_generator = tlff.ff_jones_rep_trace(
            3,
            1,
            (generator for generator in word),
        )
        self.assertExprEqual(
            trace_from_generator,
            tlff.ff_jones_rep_trace(3, 1, word),
        )

        polynomial_from_generator = tlff.ff_jones_polynomial(
            3,
            (generator for generator in word),
        )
        self.assertExprEqual(
            polynomial_from_generator,
            tlff.ff_jones_polynomial(3, word),
        )

        matrices_from_generator = tlff.ff_jones_rep_braid_word(
            3,
            1,
            (generator for generator in word),
        )
        self.assertMatricesEqual(
            matrices_from_generator,
            tlff.ff_jones_rep_braid_word(3, 1, word),
        )

        context_from_generator = tlff._ff_create_context(
            3,
            1,
            (generator for generator in word),
        )
        self.assertEqual(context_from_generator.word, tuple(word))

    def test_09b_jones_polynomial_is_stable_under_closure_moves(self):
        conjugation_cases = [
            (3, [1, -2, 1]),
            (4, [1, 2, -3, 2]),
        ]
        for n, word in conjugation_cases:
            expected = tlff.ff_jones_polynomial(n, word)
            for generator in range(1, n):
                conjugate = [generator] + word + [-generator]
                with self.subTest(move="conjugation", n=n, generator=generator):
                    actual = tlff.ff_jones_polynomial(n, conjugate)
                    self.assertExprEqual(actual, expected)

        stabilization_cases = [
            (2, [1, 1, 1]),
            (3, [1, -2, 1]),
        ]
        for n, word in stabilization_cases:
            expected = tlff.ff_jones_polynomial(n, word)
            for sign in [-1, 1]:
                stabilized_word = word + [sign * n]
                with self.subTest(move="stabilization", n=n, sign=sign):
                    actual = tlff.ff_jones_polynomial(n + 1, stabilized_word)
                    self.assertExprEqual(actual, expected)

    def test_10_ff_jones_polynomial_matches_knotinfo_random_sample(self):
        q = sp.Symbol("q")
        rows = self.load_knotinfo_fixture()
        rng = random.Random(20260413)
        sample = rng.sample(rows, 30)
        print("KnotInfo random sample:")
        for row in sample:
            print(f"  {row['name']}")

        for row in sample:
            with self.subTest(name=row["name"]):
                jones_q = sp.sympify(row["jones_q"], locals={"q": q})
                expected = sp.expand(
                    (tl.v + tl.v**-1) * jones_q.subs(q, tl.v**-2)
                )
                actual = tlff.ff_jones_polynomial(row["strands"], row["braid_word"])
                self.assertExprEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
