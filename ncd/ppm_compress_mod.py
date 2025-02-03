# Compression application using prediction by partial matching (PPM) with arithmetic coding
#
# Note that both the compressor and decompressor need to use the same PPM context modeling logic.
# The PPM algorithm can be thought of as a powerful generalization of adaptive arithmetic coding.
#
# Copyright (c) Project Nayuki
#
# https://www.nayuki.io/page/reference-arithmetic-coding
# https://github.com/nayuki/Reference-arithmetic-coding
#

import contextlib, sys, io
from complexity.ncd import arithmeticcoding, ppmmodel
import numpy as np

python3 = sys.version_info.major >= 3


# Must be at least -1 and match ppm-decompress.py. Warning: Exponential memory usage at O(257^n).
MODEL_ORDER = 3

"""
def makestream(inp):
    accum = 0
    stream = io.BytesIO()
    for i in range(len(inp)):
        if i != 0 and i % 8 == 0:
            # Write one byte to the input stream
            stream.write(bytes(accum))
            accum = 0
        # Shift accumulator left by one
        accum = accum << 1
        accum = accum + inp[i]
    # Write whatever is left in accum at the end
    stream.write(bytes(accum))
    # Change stream position to beginning
    stream.seek(0)
    return stream
"""


def makestream(a):
    """
    Accepts a numpy array and produces a string, which is then encoded into a
    byte string.
    """
    byte_stream = io.BytesIO()
    _ = byte_stream.write("".join([str(c) for c in a]).encode("ascii"))
    byte_stream.seek(0)
    return byte_stream


def compress(inp):
    outputstream = io.BytesIO()
    # From original main function, writes bits to another byte stream
    bitout = arithmeticcoding.BitOutputStream(outputstream)

    # Set up encoder and model. In this PPM model, symbol 256 represents EOF;
    # its frequency is 1 in the order -1 context but its frequency
    # is 0 in all other contexts (which have non-negative order).
    enc = arithmeticcoding.ArithmeticEncoder(32, bitout)
    model = ppmmodel.PpmModel(MODEL_ORDER, 257, 256)
    history = []

    while True:
        # Read and encode one byte
        symbol = inp.read(1)
        if len(symbol) == 0:
            break
        symbol = symbol[0] if python3 else ord(symbol)
        encode_symbol(model, history, symbol, enc)
        model.increment_contexts(history, symbol)

        if model.model_order >= 1:
            # Prepend current symbol, dropping oldest symbol if necessary
            if len(history) == model.model_order:
                history.pop()
            history.insert(0, symbol)

    encode_symbol(model, history, 256, enc)  # EOF
    enc.finish()  # Flush remaining code bits
    return outputstream


def encode_symbol(model, history, symbol, enc):
    # Try to use highest order context that exists based on the history suffix, such
    # that the next symbol has non-zero frequency. When symbol 256 is produced at a context
    # at any non-negative order, it means "escape to the next lower order with non-empty
    # context". When symbol 256 is produced at the order -1 context, it means "EOF".
    for order in reversed(range(len(history) + 1)):
        ctx = model.root_context
        for sym in history[:order]:
            assert ctx.subcontexts is not None
            ctx = ctx.subcontexts[sym]
            if ctx is None:
                break
        else:  # ctx is not None
            if symbol != 256 and ctx.frequencies.get(symbol) > 0:
                enc.write(ctx.frequencies, symbol)
                return
            # Else write context escape symbol and continue decrementing the order
            enc.write(ctx.frequencies, 256)
    # Logic for order = -1
    enc.write(model.order_minus1_freqs, symbol)


if __name__ == "__main__":
    x = np.random.binomial(1, 0.5, 2000)
    stream = makestream(x)
    output = compress(stream)
