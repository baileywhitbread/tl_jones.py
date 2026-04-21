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

- v1 used the representation $\rho\colon B_n \to \mathrm{TL}_n^\times$ defined by $\rho(\sigma_i) = v^{-1}I - v^{-2}e_i$ and the formula $J_{\widehat{\beta}}(v) = \mathrm{tr}^\mathrm{markov}(\rho(\beta))$.
- v2 uses the representation $\rho\colon B_n \to \mathrm{TL}_n^\times$ defined by $\rho(\sigma_i) = e_i - vI$ and the formula $J_{\widehat{\beta}}(v) = (-v^2)^{-w(\beta)}\,\mathrm{tr}^\mathrm{markov}(\rho(\beta))$. 
