import os

import matplotlib.pyplot as plt
import numpy

matrix_sizes = [5, 10, 25, 50, 75, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]

blas_times = [0.000041, 0.000008, 0.000082, 0.000243, 0.000458, 0.001312, 0.006275, 0.017173, 0.037515, 0.077492,
              0.127158, 0.193040, 0.279645, 0.388111, 0.553177, 0.723108, 0.921344, 1.153664, 1.414706, 1.801526,
              2.157478]


neopt_times = [0.000005, 0.000024, 0.000271, 0.002106, 0.006940, 0.016230, 0.135044, 0.470492, 1.161329, 2.275074,
               4.063601, 6.263844, 9.973713, 13.195360, 17.701904, 24.034008, 31.913750, 45.194916, 55.980133,
               75.555847, 99.743446]


opt_times = [0.000005, 0.000022, 0.000137, 0.000790, 0.002452, 0.005684, 0.042990, 0.146121, 0.353523, 0.676361,
             1.162817, 1.847915, 2.729130, 3.874825, 5.291892, 7.092333, 9.254030, 11.923586, 14.930403, 18.548702,
             23.568253]

plt.plot(matrix_sizes, blas_times, label='Blas Time')
plt.plot(matrix_sizes, neopt_times, label='Neopt Time')
plt.plot(matrix_sizes, opt_times, label='Opt Time')

plt.ylim(0, ((int) (max(neopt_times) // 10) + 2) * 10)
plt.yticks(range(0, (int) ((max(neopt_times) // 10) + 2) * 10, 5))

plt.xlabel('Matrix Size (N)')
plt.ylabel('Time (seconds)')
plt.title('Matrix Multiplication Time Performance')

plt.legend()

filename = os.path.expanduser("C:\\Users\\mblot\\Desktop\\ASC-Plots\\time_plot.png")  # set the filename to your desktop
plt.savefig(filename)

plt.show()

plt.clf()

opt_speedup = numpy.divide(neopt_times, opt_times)
blas_neopt_speedup = numpy.divide(neopt_times, blas_times)
blas_opt_speedup = numpy.divide(opt_times, blas_times)

plt.plot(matrix_sizes, opt_speedup, label='Opt Speedup over Neopt')
plt.plot(matrix_sizes, blas_neopt_speedup, label='Blas Speedup over Neopt')
plt.plot(matrix_sizes, blas_opt_speedup, label='Blas Speedup over Opt')

plt.xlabel('Matrix Size (N)')
plt.ylabel('Speedup (times)')
plt.title('Matrix Multiplication Speedup Performance')

plt.legend()

filename = os.path.expanduser("C:\\Users\\mblot\\Desktop\\ASC-Plots\\speedup_plot.png")  # set the filename to your desktop
plt.savefig(filename)

plt.show()
