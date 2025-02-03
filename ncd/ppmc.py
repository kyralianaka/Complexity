import math
import numpy as np

"""
This class creates an object from a string of ones and zeros to be compressed 
using the enclosed compress(self) method. 
The original script can be found here: https://github.com/xunzheng/ppm, but it 
has been heavily modified. The most significant changes are described below.
    - The functions were organized into a class to reinitiate the frequency table 
      for each new instance of a class (each sequence needs to be its own object).
    - The -1 order context of the table can include only '1' and '0' (not the whole
      ascii table needed for text encoding): this edit is commented and dramatically
      reduces the compression ratio but at the expense of a lower identity "score" 
      when used in the NCD formula
    - The script has been made functional in python 3.
    - The "impossible" table (whatever that is) has been initiated in the predict
      function to fix a bug.
"""


class ppmc:

    def __init__(self, arr):
        self.arr = arr
        self.bits = list()
        self.k = 3
        # If the table is a private member variable, it might change between class instances
        self.table = [dict() for i in range(self.k + 2)]  # {-1, ..., k}
        # -1 table, default ctx is ''
        self.table[0][""] = dict()
        # Instead of putting the whole ascii table in the -1 order context which would
        # be needed for text compression
        # self.table[0]['']['1'] = 1
        # self.table[0]['']['0'] = 1
        for i in range(2**8):
            self.table[0][""][chr(i)] = 1

    def predict(self, c, ctx, s, impossible=None):
        if impossible == None:
            impossible = []
        if ctx in self.table[s + 1]:
            cp = self.table[s + 1][ctx].copy()
            for key in impossible:
                if key in cp:
                    cp.pop(key)
            distinct = 0 if s == -1 else len(cp.keys())
            csum = sum(cp.values())
            if c in cp:
                prob = float(cp[c]) / float(distinct + csum)
                self.bits.append(-math.log(prob, 2))
                # print(output % ("' '" if c == ' ' else c, prob))
            else:
                if csum > 0:
                    prob = float(distinct) / float(distinct + csum)
                    self.bits.append(-math.log(prob, 2))
                    # print(output % (escape, prob))
                self.predict(c, ctx[1:], s - 1, impossible.append(cp.keys()))
        else:
            self.predict(c, ctx[1:], s - 1, impossible)

    def update(self, c, ctx):
        for i in range(0, len(ctx) + 1):
            pre = ctx[i:]
            s = len(pre)
            if pre not in self.table[s + 1]:
                self.table[s + 1][pre] = dict()
            if c not in self.table[s + 1][pre]:
                self.table[s + 1][pre][c] = 0
            self.table[s + 1][pre][c] += 1

    def compress(self):
        for i in range(len(self.arr)):
            c = self.arr[i]
            start = i - self.k if i > self.k else 0
            context = self.arr[start:i]
            self.predict(c, context, len(context), impossible=None)
            self.update(c, context)
        return sum(self.bits)


if __name__ == "__main__":
    # this is what will print if the print statements get uncommented
    escape = "<$>"
    output = "%s, %.5f"
    x = np.random.binomial(1, 0.5, 10)
    x = "".join(map(str, x))
    y = ppmc(x)
    z = y.compress()
