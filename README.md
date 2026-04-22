# tl_jones.py

Compute Jones polynomials of braid closures in Python.

Code explained in the documentation.pdf.

Example:

```python
import tl_jones as tl

# The closure of sigma_1^3 on 2 strands
print(tl.jones_polynomial(2, [1, 1, 1]))
# v**-1 + v**-3 + v**-5 - v**-9
```

#### Version history

- v1 used the braid group representation sending \sigma_i to v^{-1}I - v^{-2}e_i.
- v2 uses the braid group representation sending \sigma_i to e_i - vI.
