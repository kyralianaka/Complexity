Authors: Kyra L. Kadhim (github.com/kyralianaka), Kevin S. Brown (github.com/thelahunginjeet)

# Complexity

A python package with multiple methods of calculating the mathematical complexity of
a sequence of binary integers. These methods include Lempel-Ziv Complexity calculated
with solely python code, a sub-package for lzc calculation with a C extension that is about
50x faster, and a recently developed measure called Effort to Compress published by
Nithin Nagaraj in 2013.

In addition to these measures of a sequence's complexity, some relative measures of
complexity have been included. The minimum hamming distance between two or more sequences
can be calculating using Hamming.py, and the NCD sub-package can be used to calculate the
normalized compression distance between two sequences, or the normalized compression
distance matrix of multiple sequences.

To calculate the normalized compression distance of sequences, different compressors
may be used. The compressors supported by NCD.py include gzip, bzip, snappy, ppm
(prediction by partial matching algorithm developed by [Nayuki](https://github.com/nayuki/Reference-arithmetic-coding)), ppmc, and the LZ Complexity itself substituted for
the compression value.

## Installation

This package may be installed on both Windows and UNIX operating systems provided they are
using python 3. The pip package manager is recommended for installation.
```bash
pip install complexity
```
**Note:** Your machine may already have a package installed named complexity. If so,
it will be necessary to navigate to the complexity directory (after cloning the repository)
and use the following.
```bash
pip install .
```

## Usage

To calculate the Lempel-Ziv Complexity of a sequence using the C extension:
```python
from complexity.lzc.lzc import lz_complexity
lz_complexity([1,0,1,0,0]) # returns 3
```

To calculate the Lempel-Ziv Complexity of a sequence using python:
```python
from complexity.lz_complexity_python import lz_complexity
lz_complexity([1,0,1,0,0,1,1]) # returns 3
```

To calculate the Effort to Compress value of a sequence:
```python
from complexity.etc import etc
etc([0,1,0,0,1,0,1,1]) # returns ([8.721, 8.297, 7.443, 6.198, 2.721, 0.0], 5)
```

To calculate the NCD between sequences:
```python
import numpy as np
import gzip
from complexity.ncd import NCD
x = np.tile([1,0,1,0],10)
y = np.tile([1,1,0,1],10)
NCD.NCD_pairwise(x,y,gzip) # returns 0.206
```

## License

See license text file included in the package.
