import numpy as np
import matplotlib.pyplot as plt

N = 61
strain = np.zeros([N - 1, 6])
T = np.zeros(N - 1)

for i in range(1, N):
    eta_h = np.fromfile(f'{i}/monitor_eta_h', dtype=float)
    eta_h = np.reshape(eta_h, (-1, 6))
    eta_h_avg = np.average(eta_h, axis=0)
    strain[i - 1, :] = eta_h_avg * 1e2
    T[i - 1] = 350 - 5 * (i - 1)

T = T[::-1]
strain = strain[::-1, :]
strain[:, 4] *= -1
strain[:, 5] *= -1

idx_strain = [0, 2, 1, 3, 5, 4]
for i in range(6):
    plt.plot(T, strain[:, idx_strain[i]], label=f'$\\eta_{i + 1}$')
plt.xlabel('Temperature (K)')
plt.ylabel('Strain ($10^{-2}$)')
plt.legend(loc=(1/7,1/4))
plt.tight_layout()
plt.savefig('fig.png', dpi=300)

