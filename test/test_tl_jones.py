import ast
import csv
import random
import sys
import unittest
from math import comb
from pathlib import Path

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import tl_jones as tl

FIXTURE_PATH = PROJECT_ROOT / "data" / "jones_braid_knots.csv"


class JonesPolynomialTests(unittest.TestCase):
    def test_01_validation_errors(self):
        self.assertEqual(tl.link_states(0, 0), ((),))
        with self.assertRaises(ValueError):
            tl.link_states(3, 2)
        with self.assertRaises(ValueError):
            tl.jones_rep_braid_generator(3, 0, 0)
        with self.assertRaises(ValueError):
            tl.jones_rep_braid_generator(3, 0, 3)
        with self.assertRaises(ValueError):
            tl.jones_polynomial(2, [0])
        with self.assertRaises(ValueError):
            tl.jones_polynomial(2, [2])
        with self.assertRaises(ValueError):
            tl.jones_polynomial(2, [1.0])

    def test_01a_gnf_to_braid_word_examples(self):
        self.assertEqual(tl.gnf_to_braid_word(1, 1, ()), [])
        self.assertEqual(tl.gnf_to_braid_word(3, 1, ()), [1, 2, 1])
        self.assertEqual(tl.gnf_to_braid_word(4, 1, ()), [1, 2, 3, 1, 2, 1])
        self.assertEqual(tl.gnf_to_braid_word(3, -1, ()), [-1, -2, -1])
        self.assertEqual(tl.gnf_to_braid_word(3, 0, ((2, 1, 3),)), [1])
        self.assertEqual(tl.gnf_to_braid_word(3, 0, ((2, 3, 1),)), [1, 2])
        self.assertEqual(tl.gnf_to_braid_word(3, 0, ((3, 1, 2),)), [2, 1])
        self.assertEqual(
            tl.gnf_to_braid_word(3, -1, ((2, 1, 3),)),
            [-1, -2, -1, 1],
        )
        self.assertEqual(
            tl.gnf_to_braid_word(3, 0, ((2, 1, 3), (2, 3, 1))),
            [1, 1, 2],
        )

    def test_01b_gnf_to_braid_word_validation_errors(self):
        with self.assertRaises(ValueError):
            tl.gnf_to_braid_word(0, 0, ())
        with self.assertRaises(ValueError):
            tl.gnf_to_braid_word(3, 1.0, ())
        with self.assertRaises(ValueError):
            tl.gnf_to_braid_word(3, 0, ((2, 1),))
        with self.assertRaises(ValueError):
            tl.gnf_to_braid_word(3, 0, ((2, 2, 3),))
        with self.assertRaises(ValueError):
            tl.gnf_to_braid_word(3, 0, ((1, 2, 3),))
        with self.assertRaises(ValueError):
            tl.gnf_to_braid_word(3, 0, ((3, 2, 1),))
        with self.assertRaises(ValueError):
            tl.gnf_to_braid_word(3, 0, ((2, 1, 3), (1, 3, 2)))

    def test_02_projlen_hard_coded_cases(self):
        cases = [
            (tl.v + tl.v**-1, 2),
            (tl.v**-100 + tl.v, 101),
            (tl.v**17, 0),
            (5, 0),
            (0, 0),
            (3 * tl.v**4 - 2 * tl.v**-3 + 7, 7),
            (tl.v**8 - tl.v**3 + 2, 8),
            (-4 * tl.v**-7 + tl.v**-2, 5),
            (2 * tl.v**9 + 3 * tl.v**-9, 18),
            (tl.v**5 + tl.v**4 + tl.v**3, 2),
            (-tl.v**12 + 6 * tl.v**12, 0),
            (7 * tl.v**-1 - 3 * tl.v**-6 + 2 * tl.v**2, 8),
            (tl.v**20 - tl.v**19 + tl.v**-20, 40),
            (-2 * tl.v**15 + 9 * tl.v**7 - tl.v**6 + 1, 15),
            (tl.v**-11 - 4 * tl.v**-13 - 2 * tl.v**-12, 2),
            (3 * tl.v**14 - 5 * tl.v**-4 + 8 * tl.v**2 - 1, 18),
            (tl.v - tl.v, 0),
            (2 * tl.v**-25 + tl.v**25, 50),
            (-3 * tl.v**10 + tl.v**8 - tl.v**2, 8),
            (6 * tl.v**-3 + 4 * tl.v**-3 + 1, 3),
        ]

        for polynomial, expected in cases:
            with self.subTest(polynomial=polynomial):
                self.assertEqual(tl.projlen(polynomial), expected)

    def test_03_link_states_examples(self):
        examples = [
            (1, 0, ((-1,),)),
            (2, 1, ((1, 0),)),
            (3, 1, ((-1, 2, 1), (1, 0, -1))),
            (4, 1, ((-1, -1, 3, 2), (-1, 2, 1, -1), (1, 0, -1, -1))),
            (4, 2, ((1, 0, 3, 2), (3, 2, 1, 0))),
            (
                5,
                2,
                (
                    (-1, 2, 1, 4, 3),
                    (-1, 4, 3, 2, 1),
                    (1, 0, -1, 4, 3),
                    (1, 0, 3, 2, -1),
                    (3, 2, 1, 0, -1),
                ),
            ),
        ]

        for n, r, expected in examples:
            with self.subTest(n=n, r=r):
                self.assertEqual(tl.link_states(n, r), expected)

    def test_04_link_state_dimensions(self):
        for n in range(1, 11):
            for r in range(n // 2 + 1):
                previous = 0 if r == 0 else comb(n, r - 1)
                expected = comb(n, r) - previous
                self.assertEqual(len(tl.link_states(n, r)), expected)

    def test_05_temperley_lieb_relations(self):
        for n in range(2, 11):
            for r in range(n // 2 + 1):
                generators = [
                    tl._temperley_lieb_matrix(n, r, i) for i in range(1, n)
                ]

                for index, generator in enumerate(generators):
                    left = generator * generator
                    right = tl.quantum_two * generator
                    self.assertEqual(left.shape, right.shape)
                    for row in range(left.rows):
                        for column in range(left.cols):
                            self.assertEqual(
                                sp.simplify(left[row, column] - right[row, column]),
                                0,
                            )

                    for other_index, other_generator in enumerate(generators):
                        if abs(index - other_index) > 1:
                            left = generator * other_generator
                            right = other_generator * generator
                            self.assertEqual(left.shape, right.shape)
                            for row in range(left.rows):
                                for column in range(left.cols):
                                    self.assertEqual(
                                        sp.simplify(
                                            left[row, column] - right[row, column]
                                        ),
                                        0,
                                    )

                    if index + 1 < len(generators):
                        left = generator * generators[index + 1] * generator
                        right = generator
                        self.assertEqual(left.shape, right.shape)
                        for row in range(left.rows):
                            for column in range(left.cols):
                                self.assertEqual(
                                    sp.simplify(left[row, column] - right[row, column]),
                                    0,
                                )

    def test_05a_braid_generator_formula(self):
        for n in range(2, 6):
            for r in range(n // 2 + 1):
                dimension = len(tl.link_states(n, r))
                identity = sp.eye(dimension)
                for i in range(1, n):
                    positive_expected = (
                        tl._temperley_lieb_matrix(n, r, i) - tl.v * identity
                    )
                    positive_actual = tl.jones_rep_braid_generator(n, r, i)
                    self.assertEqual(positive_actual, positive_expected)

                    inverse_expected = (
                        tl._temperley_lieb_matrix(n, r, i) - tl.v**-1 * identity
                    )
                    inverse_actual = tl.jones_rep_braid_generator(
                        n, r, i, inverse=True
                    )
                    self.assertEqual(inverse_actual, inverse_expected)

    def test_06_braid_inverse_relation(self):
        for n in range(2, 11):
            for r in range(n // 2 + 1):
                identity = sp.eye(len(tl.link_states(n, r)))
                for i in range(1, n):
                    positive = tl.jones_rep_braid_generator(n, r, i)
                    negative = tl.jones_rep_braid_generator(
                        n, r, i, inverse=True
                    )
                    self.assertEqual(negative, tl.jones_rep_braid_word(n, r, [-i]))
                    left = positive * negative
                    right = identity
                    self.assertEqual(left.shape, right.shape)
                    for row in range(left.rows):
                        for column in range(left.cols):
                            self.assertEqual(
                                sp.simplify(left[row, column] - right[row, column]),
                                0,
                            )

                    left = negative * positive
                    right = identity
                    self.assertEqual(left.shape, right.shape)
                    for row in range(left.rows):
                        for column in range(left.cols):
                            self.assertEqual(
                                sp.simplify(left[row, column] - right[row, column]),
                                0,
                            )

    def test_07_braid_commuting_relation(self):
        for n in range(4, 11):
            for r in range(n // 2 + 1):
                for i in range(1, n):
                    for j in range(1, n):
                        if abs(i - j) > 1:
                            left = (
                                tl.jones_rep_braid_generator(n, r, i)
                                * tl.jones_rep_braid_generator(n, r, j)
                            )
                            right = (
                                tl.jones_rep_braid_generator(n, r, j)
                                * tl.jones_rep_braid_generator(n, r, i)
                            )
                            self.assertEqual(left.shape, right.shape)
                            for row in range(left.rows):
                                for column in range(left.cols):
                                    self.assertEqual(
                                        sp.simplify(
                                            left[row, column] - right[row, column]
                                        ),
                                        0,
                                    )

    def test_08_braid_relation(self):
        for n in range(3, 11):
            for r in range(n // 2 + 1):
                for i in range(1, n - 1):
                    sigma_i = tl.jones_rep_braid_generator(n, r, i)
                    sigma_next = tl.jones_rep_braid_generator(n, r, i + 1)
                    left = sigma_i * sigma_next * sigma_i
                    right = sigma_next * sigma_i * sigma_next
                    self.assertEqual(left.shape, right.shape)
                    for row in range(left.rows):
                        for column in range(left.cols):
                            self.assertEqual(
                                sp.simplify(left[row, column] - right[row, column]),
                                0,
                            )

    def test_09_jones_polynomial_simple_cases(self):
        self.assertEqual(
            sp.simplify(tl.jones_polynomial(1, []) - (tl.v + tl.v**-1)),
            0,
        )

        for n in range(1, 11):
            self.assertEqual(
                sp.simplify(
                    tl.jones_polynomial(n, []) - sp.expand((tl.v + tl.v**-1) ** n)
                ),
                0,
            )

        self.assertEqual(
            sp.simplify(tl.jones_polynomial(2, [1]) - (tl.v + tl.v**-1)),
            0,
        )
        self.assertEqual(
            sp.simplify(tl.jones_polynomial(2, [-1]) - (tl.v + tl.v**-1)),
            0,
        )
        self.assertEqual(
            sp.simplify(
                tl.jones_polynomial(2, [1, 1])
                - (1 + tl.v**-2 + tl.v**-4 + tl.v**-6)
            ),
            0,
        )
        self.assertEqual(
            sp.simplify(
                tl.jones_polynomial(2, [1, 1, 1])
                - (tl.v**-1 + tl.v**-3 + tl.v**-5 - tl.v**-9)
            ),
            0,
        )

    def test_10_jones_polynomial_matches_knotinfo_random_sample(self):
        # The fixture stores normalized Jones polynomials in q, while this
        # module returns the unnormalized Markov trace in v; the matching
        # convention is J_ours(v) = (v + v^-1) V_KnotInfo(q) at q = v^-2.
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
                actual = tl.jones_polynomial(row["strands"], row["braid_word"])
                self.assertEqual(sp.simplify(actual - expected), 0)

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


if __name__ == "__main__":
    unittest.main()
