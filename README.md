# tl_jones.py

Compute Jones polynomials of braid closures in Python.

Code explained in the documentation.pdf.

Example:

```python
import tl_jones as tl

# The closure of sigma_1^3 on 2 strands
print(tl.jones_polynomial(2, [1, 1, 1]))
# v**-1 + v**-3 + v**-5 - v**-9
