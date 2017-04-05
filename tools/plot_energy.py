#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
import sys

class Energy:
    def __init__(self, fname):
        self.fromfile(fname)

    def fromfile(self, fname):
        f = open(fname, 'rb')

        self.magic         = f.read(3)
        self.version       = self.read_int(f)
        self.revision      = self.read_int(f)
        self.endianness    = self.read_int(f)
        self.header_length = self.read_int(f)
        self.num_sz        = self.read_int(f)
        self.en_nvars      = self.read_int(f)
        self.id_length     = self.read_int(f)
        self.varnames      = f.read(self.en_nvars*self.id_length).split()

        self.data = np.fromfile(f, dtype=np.float64)
        n_timesteps = len(self.data)/self.en_nvars
        self.data = self.data.reshape((n_timesteps, self.en_nvars))

        f.close()

    def plot(self, index):
        plt.semilogy(en.data[:, 0], en.data[:,index], label=str(index))

    def fortran_to_int(self, n):
        return int(n[::-1].encode('hex'), 16)

    def read_int(self, fp):
        integer_byte_size = 4
        fortran_int = fp.read(integer_byte_size)
        return self.fortran_to_int(fortran_int)

fname = sys.argv[1]
en = Energy(fname)
print(en.varnames)
fig, ax = plt.subplots()
en.plot(4)
if en.en_nvars > 6:
    en.plot(8)
    en.plot(9)
ax.set_ylim(1e-12, 1e-2)
ax.set_xlim(1, 10)
plt.show()
