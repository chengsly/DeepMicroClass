import numpy as np
import sys


arr = np.load(sys.argv[1], allow_pickle=True)
print("shape: {}".format(arr.shape))
print("size: {}".format(float(arr.nbytes/1000000000)))

