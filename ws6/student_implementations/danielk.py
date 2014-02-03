import numpy as np

def dft(x):
  '''
  Calculates the discrete Fourier transform of x.
  '''
  n = len(x)
  w = np.zeros((n, n),dtype=complex)
  for j in range(n):
    for k in range(n):
      w[j, k] = np.exp(-2 * np.pi * 1j * j * k / n)
  return w.dot(x)
