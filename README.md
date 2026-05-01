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

- v1.1 swapped the braid group representation; instead of sending \sigma_i to v^{-1}I - v^{-2}e_i, we send \sigma_i to e_i - vI.
- v2.0 added tl_jones_ff.py which implements inverse discrete Fourier transform to speed up calculations.
