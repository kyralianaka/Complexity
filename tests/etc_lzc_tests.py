from complexity.etc import etc
from complexity.lzc.lzc import lz_complexity
import numpy as np
import matplotlib.pyplot as plt

# Plot to show lzc and etc dependence on length for random sequences and
# plateu of both lzc and etc for a repeating sequence at the length that
# the sequence starts to repeat itself
x = np.random.binomial(1, 0.5, 1000)  # random seq
y = np.random.binomial(1, 0.5, 500)
y = np.tile(y, 2)  # repeated seq

lens = np.arange(100, 1100, step=100)
etc_x = np.zeros_like(lens)
etc_y = np.zeros_like(lens)
lzc_x = np.zeros_like(lens)
lzc_y = np.zeros_like(lens)

for idx, l in enumerate(lens):
    etc_x[idx] = etc(x[:l])[1]
    etc_y[idx] = etc(y[:l])[1]
    lzc_x[idx] = lz_complexity(x[:l])
    lzc_y[idx] = lz_complexity(y[:l])

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].plot(lens, etc_x, label="ETC", color="b")
axs[0].plot(lens, lzc_x, label="LZC", color="r")
axs[1].plot(lens, etc_y, label="ETC", color="b")
axs[1].plot(lens, lzc_y, label="LZC", color="r")

axs[0].title.set_text("Random Sequence")
axs[1].title.set_text("Repeated Sequence (@n=500)")
plt.setp(axs, xlabel="Sequence Length", ylabel="Complexity")

plt.legend()
plt.show()
