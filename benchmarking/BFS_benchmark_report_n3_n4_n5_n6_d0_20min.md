# BFS benchmark report for `d = 0` in `B_3`, `B_4`, `B_5`, and `B_6`

Generated: `2026-05-01T08:58:48+10:00`

This report compares the symbolic implementation `tl_jones.py` with the finite-field implementation `tl_jones_ff.py` on breadth-first Garside-normal-form trees.

## Setup

- Time budget: `20` minutes per braid group per code.
- Protocol: each cell was run in a fresh Python process, sequentially.
- Timed section: starts after BFS tree metadata is built; Python caches are cold for each cell.
- Workload at each node: build `gnf_to_braid_word(n, 0, factors)` and compute the Jones polynomial.
- Depth convention: the empty braid is depth `0`; depth equals the number of simple factors.
- Tree rule: simple factors exclude identity and `Delta`; transitions require `LeftDescentSet(next) subseteq RightDescentSet(previous)`.
- Timeout accounting: timeout is checked after each completed braid, so elapsed time can exceed the nominal budget by one computation.
- Raw JSON files: `benchmark_results`.
- Python: `3.12.0 (tags/v3.12.0:0fb18b0, Oct  2 2023, 13:03:39) [MSC v.1935 64 bit (AMD64)]`.
- SymPy: `1.13.1`.
- Platform: `Windows-11-10.0.26200-SP0`.

## Results

### Processed braids by depth layer

![Processed braids by BFS depth layer](benchmark_figures/bfs_processed_braids_by_depth_layer.svg)


### Completed depth before timeout

| Implementation | `B_3` | `B_4` | `B_5` | `B_6` |
| --- | ---: | ---: | ---: | ---: |
| `tl_jones.py` | 13 | 5 | 2 | 1 |
| `tl_jones_ff.py` | 15 | 5 | 3 | 1 |

### Throughput summary

| Implementation | Group | Complete depth | Timeout depth | Total processed braids | Complete-depth braids | Overall braids/s | Timeout-layer braids/s | Timeout-layer completion |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `tl_jones.py` | `B_3` | 13 | 14 | 38,320 | 32,765 | 31.93 | 19.19 | 16.95% |
| `tl_jones_ff.py` | `B_3` | 15 | 16 | 158,665 | 131,069 | 132.22 | 99.48 | 21.05% |
| `tl_jones.py` | `B_4` | 5 | 6 | 54,947 | 37,175 | 45.79 | 28.05 | 10.65% |
| `tl_jones_ff.py` | `B_4` | 5 | 6 | 137,841 | 37,175 | 114.87 | 102.25 | 60.35% |
| `tl_jones.py` | `B_5` | 2 | 3 | 33,799 | 3,531 | 28.16 | 25.80 | 42.06% |
| `tl_jones_ff.py` | `B_5` | 3 | 4 | 91,282 | 75,489 | 76.07 | 73.51 | 1.13% |
| `tl_jones.py` | `B_6` | 1 | 2 | 23,529 | 719 | 19.61 | 19.09 | 25.49% |
| `tl_jones_ff.py` | `B_6` | 1 | 2 | 46,549 | 719 | 38.79 | 38.36 | 51.22% |


### Nodes per depth


| Depth | `B_3` | `B_4` | `B_5` | `B_6` |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | 1 | 1 | 1 |
| 1 | 4 | 22 | 118 | 718 |
| 2 | 8 | 164 | 3,412 | 89,482 |
| 3 | 16 | 982 | 71,958 | n/a |
| 4 | 32 | 5,528 | 1,394,072 | n/a |
| 5 | 64 | 30,478 | n/a | n/a |
| 6 | 128 | 166,796 | n/a | n/a |
| 7 | 256 | n/a | n/a | n/a |
| 8 | 512 | n/a | n/a | n/a |
| 9 | 1,024 | n/a | n/a | n/a |
| 10 | 2,048 | n/a | n/a | n/a |
| 11 | 4,096 | n/a | n/a | n/a |
| 12 | 8,192 | n/a | n/a | n/a |
| 13 | 16,384 | n/a | n/a | n/a |
| 14 | 32,768 | n/a | n/a | n/a |
| 15 | 65,536 | n/a | n/a | n/a |
| 16 | 131,072 | n/a | n/a | n/a |


## Dense summary table

### `B_3`

| Implementation | Result |
| --- | --- |
| `tl_jones.py` | Depth 13<br>Completed time: 910.615 s<br>Braids through depth: 32,765<br>Final progress: depth 14: 5,555 / 32,768 (16.953%)<br>Final-depth time: 289.414 s<br>Final-depth braids: 5,555 |
| `tl_jones_ff.py` | Depth 15<br>Completed time: 922.594 s<br>Braids through depth: 131,069<br>Final progress: depth 16: 27,596 / 131,072 (21.054%)<br>Final-depth time: 277.409 s<br>Final-depth braids: 27,596 |

### `B_4`

| Implementation | Result |
| --- | --- |
| `tl_jones.py` | Depth 5<br>Completed time: 566.468 s<br>Braids through depth: 37,175<br>Final progress: depth 6: 17,772 / 166,796 (10.655%)<br>Final-depth time: 633.531 s<br>Final-depth braids: 17,772 |
| `tl_jones_ff.py` | Depth 5<br>Completed time: 215.515 s<br>Braids through depth: 37,175<br>Final progress: depth 6: 100,666 / 166,796 (60.353%)<br>Final-depth time: 984.493 s<br>Final-depth braids: 100,666 |

### `B_5`

| Implementation | Result |
| --- | --- |
| `tl_jones.py` | Depth 2<br>Completed time: 26.715 s<br>Braids through depth: 3,531<br>Final progress: depth 3: 30,268 / 71,958 (42.063%)<br>Final-depth time: 1,173.346 s<br>Final-depth braids: 30,268 |
| `tl_jones_ff.py` | Depth 3<br>Completed time: 985.180 s<br>Braids through depth: 75,489<br>Final progress: depth 4: 15,793 / 1,394,072 (1.133%)<br>Final-depth time: 214.827 s<br>Final-depth braids: 15,793 |

### `B_6`

| Implementation | Result |
| --- | --- |
| `tl_jones.py` | Depth 1<br>Completed time: 4.906 s<br>Braids through depth: 719<br>Final progress: depth 2: 22,810 / 89,482 (25.491%)<br>Final-depth time: 1,195.090 s<br>Final-depth braids: 22,810 |
| `tl_jones_ff.py` | Depth 1<br>Completed time: 5.184 s<br>Braids through depth: 719<br>Final progress: depth 2: 45,830 / 89,482 (51.217%)<br>Final-depth time: 1,194.812 s<br>Final-depth braids: 45,830 |

## Detailed timings

### `B_3`

#### `tl_jones.py`

- Backend key: `symbolic`
- Timed out: `True`
- Completed depth: `13`
- Total processed nodes: `38,320`
- Simple factors: `4`
- Non-root child counts: min `2`, max `2`, avg `2.000`
- Setup time excluded from timed section: `0.000 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.002 s | 0.003 s | complete |
| 1 | 4 | 4 | 0.006 s | 0.011 s | complete |
| 2 | 8 | 8 | 0.011 s | 0.023 s | complete |
| 3 | 16 | 16 | 0.029 s | 0.055 s | complete |
| 4 | 32 | 32 | 0.056 s | 0.113 s | complete |
| 5 | 64 | 64 | 0.155 s | 0.271 s | complete |
| 6 | 128 | 128 | 0.456 s | 0.729 s | complete |
| 7 | 256 | 256 | 1.078 s | 1.810 s | complete |
| 8 | 512 | 512 | 2.628 s | 4.441 s | complete |
| 9 | 1,024 | 1,024 | 6.097 s | 10.540 s | complete |
| 10 | 2,048 | 2,048 | 15.157 s | 25.699 s | complete |
| 11 | 4,096 | 4,096 | 50.779 s | 76.481 s | complete |
| 12 | 8,192 | 8,192 | 163.135 s | 239.619 s | complete |
| 13 | 16,384 | 16,384 | 670.994 s | 910.615 s | complete |
| 14 | 32,768 | 5,555 | 289.414 s | 1,200.032 s | partial |

#### `tl_jones_ff.py`

- Backend key: `ff`
- Timed out: `True`
- Completed depth: `15`
- Total processed nodes: `158,665`
- Simple factors: `4`
- Non-root child counts: min `2`, max `2`, avg `2.000`
- Setup time excluded from timed section: `0.000 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.004 s | 0.006 s | complete |
| 1 | 4 | 4 | 0.004 s | 0.012 s | complete |
| 2 | 8 | 8 | 0.007 s | 0.021 s | complete |
| 3 | 16 | 16 | 0.013 s | 0.035 s | complete |
| 4 | 32 | 32 | 0.027 s | 0.065 s | complete |
| 5 | 64 | 64 | 0.067 s | 0.134 s | complete |
| 6 | 128 | 128 | 0.173 s | 0.310 s | complete |
| 7 | 256 | 256 | 0.448 s | 0.761 s | complete |
| 8 | 512 | 512 | 1.162 s | 1.925 s | complete |
| 9 | 1,024 | 1,024 | 2.731 s | 4.659 s | complete |
| 10 | 2,048 | 2,048 | 6.262 s | 10.923 s | complete |
| 11 | 4,096 | 4,096 | 16.227 s | 27.153 s | complete |
| 12 | 8,192 | 8,192 | 36.053 s | 63.208 s | complete |
| 13 | 16,384 | 16,384 | 87.355 s | 150.566 s | complete |
| 14 | 32,768 | 32,768 | 216.111 s | 366.680 s | complete |
| 15 | 65,536 | 65,536 | 555.912 s | 922.594 s | complete |
| 16 | 131,072 | 27,596 | 277.409 s | 1,200.007 s | partial |

### `B_4`

#### `tl_jones.py`

- Backend key: `symbolic`
- Timed out: `True`
- Completed depth: `5`
- Total processed nodes: `54,947`
- Simple factors: `22`
- Non-root child counts: min `3`, max `11`, avg `7.455`
- Setup time excluded from timed section: `0.000 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.002 s | 0.004 s | complete |
| 1 | 22 | 22 | 0.042 s | 0.047 s | complete |
| 2 | 164 | 164 | 0.445 s | 0.495 s | complete |
| 3 | 982 | 982 | 5.658 s | 6.156 s | complete |
| 4 | 5,528 | 5,528 | 53.126 s | 59.286 s | complete |
| 5 | 30,478 | 30,478 | 507.179 s | 566.468 s | complete |
| 6 | 166,796 | 17,772 | 633.531 s | 1,200.002 s | partial |

#### `tl_jones_ff.py`

- Backend key: `ff`
- Timed out: `True`
- Completed depth: `5`
- Total processed nodes: `137,841`
- Simple factors: `22`
- Non-root child counts: min `3`, max `11`, avg `7.455`
- Setup time excluded from timed section: `0.000 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.002 s | 0.005 s | complete |
| 1 | 22 | 22 | 0.022 s | 0.029 s | complete |
| 2 | 164 | 164 | 0.207 s | 0.238 s | complete |
| 3 | 982 | 982 | 2.181 s | 2.421 s | complete |
| 4 | 5,528 | 5,528 | 20.354 s | 22.778 s | complete |
| 5 | 30,478 | 30,478 | 192.734 s | 215.515 s | complete |
| 6 | 166,796 | 100,666 | 984.493 s | 1,200.012 s | partial |

### `B_5`

#### `tl_jones.py`

- Backend key: `symbolic`
- Timed out: `True`
- Completed depth: `2`
- Total processed nodes: `33,799`
- Simple factors: `118`
- Non-root child counts: min `4`, max `59`, avg `28.915`
- Setup time excluded from timed section: `0.001 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.003 s | 0.006 s | complete |
| 1 | 118 | 118 | 0.307 s | 0.314 s | complete |
| 2 | 3,412 | 3,412 | 26.399 s | 26.715 s | complete |
| 3 | 71,958 | 30,268 | 1,173.346 s | 1,200.064 s | partial |

#### `tl_jones_ff.py`

- Backend key: `ff`
- Timed out: `True`
- Completed depth: `3`
- Total processed nodes: `91,282`
- Simple factors: `118`
- Non-root child counts: min `4`, max `59`, avg `28.915`
- Setup time excluded from timed section: `0.001 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.003 s | 0.006 s | complete |
| 1 | 118 | 118 | 0.199 s | 0.207 s | complete |
| 2 | 3,412 | 3,412 | 17.975 s | 18.185 s | complete |
| 3 | 71,958 | 71,958 | 966.993 s | 985.180 s | complete |
| 4 | 1,394,072 | 15,793 | 214.827 s | 1,200.010 s | partial |

### `B_6`

#### `tl_jones.py`

- Backend key: `symbolic`
- Timed out: `True`
- Completed depth: `1`
- Total processed nodes: `23,529`
- Simple factors: `718`
- Non-root child counts: min `5`, max `359`, avg `124.627`
- Setup time excluded from timed section: `0.033 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.011 s | 0.014 s | complete |
| 1 | 718 | 718 | 4.884 s | 4.906 s | complete |
| 2 | 89,482 | 22,810 | 1,195.090 s | 1,200.002 s | partial |

#### `tl_jones_ff.py`

- Backend key: `ff`
- Timed out: `True`
- Completed depth: `1`
- Total processed nodes: `46,549`
- Simple factors: `718`
- Non-root child counts: min `5`, max `359`, avg `124.627`
- Setup time excluded from timed section: `0.032 s`

| Depth | Nodes in layer | Processed | Level time | Cumulative time | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 1 | 0.004 s | 0.006 s | complete |
| 1 | 718 | 718 | 5.175 s | 5.184 s | complete |
| 2 | 89,482 | 45,830 | 1,194.812 s | 1,200.001 s | partial |
