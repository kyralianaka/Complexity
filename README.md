Authors: Kyra L. Kadhim (github.com/kyralianaka), Kevin S. Brown (github.com/thelahunginjeet)

# Complexity

A Python package with multiple methods of calculating the mathematical complexity of
a sequence of binary integers. These methods include Lempel-Ziv Complexity calculated
with solely python code, a sub-package for LZC calculation with a C extension that is about
50x faster, and a recently developed measure called Effort to Compress published by
Nithin Nagaraj in 2013.

In addition to these measures of a sequence's complexity, some relative measures of
complexity have been included. The minimum hamming distance between two or more sequences
can be calculating using hamming.py, and the NCD sub-package can be used to calculate the
Normalized Compression Distance between two sequences, or the Normalized Compression
Distance matrix of multiple sequences.

To calculate the Normalized Compression Distance of sequences, different compressors
may be used. The compressors supported by NCD.py include gzip, bzip, snappy, ppm
(prediction by partial matching algorithm developed by [Nayuki](https://github.com/nayuki/Reference-arithmetic-coding)), ppmc, and the LZ Complexity itself substituted for
the compressed length.

## Installation

This package may be installed on both Windows and UNIX operating systems provided they are
using Python 3. The pip package manager is recommended for installation.
```bash
pip install complexity
```
**Note:** There are other packages named complexity (sorry), so please configure your environment accordingly.

## Usage

To calculate the Lempel-Ziv Complexity of a sequence using the C extension:
```python
from complexity.lzc.lzc import lz_complexity
lz_complexity([1, 0, 1, 0, 0]) # returns 3
```

To calculate the Lempel-Ziv Complexity of a sequence using Python:
```python
from complexity.lz_complexity_python import lz_complexity
lz_complexity([1, 0, 1, 0, 0, 1, 1]) # returns 3
```

To calculate the Effort to Compress value of a sequence:
```python
from complexity.etc import etc
etc([0, 1, 0, 0, 1, 0, 1, 1]) # returns ([8.721, 8.297, 7.443, 6.198, 2.721, 0.0], 5)
```

To calculate the NCD between sequences:
```python
import numpy as np
import gzip
from complexity.ncd import NCD
x = np.tile([1, 0, 1, 0], 10)
y = np.tile([1, 1, 0, 1], 10)
NCD.NCD_pairwise(x, y, gzip) # returns 0.206
```

## License

See license text file included in the package.
