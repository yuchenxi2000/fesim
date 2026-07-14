"""
Python reader of dipole interaction matrix
"""
import struct
import numpy as np


def array_index(i, j, k, Nx, Ny, Nz):
    return Nz * (Ny * i + j) + k


class QMat:
    def __init__(self, filepath):
        with open(filepath, 'rb') as fin:
            file_version, = struct.unpack('i', fin.read(4))
            if file_version != 2:
                raise Exception(f'unsupported file version {file_version}!')
            header_length, = struct.unpack('i', fin.read(4))
            self.Nx, self.Ny, self.Nz, _align_4, \
                self.ax, self.ay, self.az, \
                self.lambda_, self.q_tol = struct.unpack('iiiiddddd', fin.read(header_length))
            self.Nx_h = (self.Nx + 2) // 2
            self.Ny_h = (self.Ny + 2) // 2
            self.Nz_h = (self.Nz + 2) // 2
            self.N_h = self.Nx_h * self.Ny_h * self.Nz_h
            self.qxx = np.frombuffer(fin.read(self.N_h * 8), dtype=float)
            self.qyy = np.frombuffer(fin.read(self.N_h * 8), dtype=float)
            self.qzz = np.frombuffer(fin.read(self.N_h * 8), dtype=float)
            self.qyz = np.frombuffer(fin.read(self.N_h * 8), dtype=float)
            self.qzx = np.frombuffer(fin.read(self.N_h * 8), dtype=float)
            self.qxy = np.frombuffer(fin.read(self.N_h * 8), dtype=float)
            self.Qxx = self.convert_q_to_Q0(self.qxx, [1.0, 1.0, 1.0])
            self.Qyy = self.convert_q_to_Q0(self.qyy, [1.0, 1.0, 1.0])
            self.Qzz = self.convert_q_to_Q0(self.qzz, [1.0, 1.0, 1.0])
            self.Qyz = self.convert_q_to_Q0(self.qyz, [1.0, -1.0, -1.0])
            self.Qzx = self.convert_q_to_Q0(self.qzx, [-1.0, 1.0, -1.0])
            self.Qxy = self.convert_q_to_Q0(self.qxy, [-1.0, -1.0, 1.0])
            # Precompute FFT of each Q component for fast convolution
            self._Qxx_f = np.fft.fftn(self.Qxx)
            self._Qyy_f = np.fft.fftn(self.Qyy)
            self._Qzz_f = np.fft.fftn(self.Qzz)
            self._Qyz_f = np.fft.fftn(self.Qyz)
            self._Qzx_f = np.fft.fftn(self.Qzx)
            self._Qxy_f = np.fft.fftn(self.Qxy)

    def convert_q_to_Q0(self, q, sign_change):
        Q0 = np.zeros([self.Nx, self.Ny, self.Nz], dtype=float)
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    sign = 1.0
                    if i < self.Nx_h:
                        i2 = i
                    else:
                        i2 = self.Nx - i
                        sign *= sign_change[0]
                    if j < self.Ny_h:
                        j2 = j
                    else:
                        j2 = self.Ny - j
                        sign *= sign_change[1]
                    if k < self.Nz_h:
                        k2 = k
                    else:
                        k2 = self.Nz - k
                        sign *= sign_change[2]
                    Q0[i, j, k] = q[array_index(i2, j2, k2, self.Nx_h, self.Ny_h, self.Nz_h)] * sign
        return Q0

    def dp_int_direct(self, dp: np.ndarray, normalize: bool = False) -> float:
        """Translated from C++ code, test only"""
        s = 0.0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    tmp = np.zeros(3, dtype=float)
                    for i1 in range(self.Nx):
                        for j1 in range(self.Ny):
                            for k1 in range(self.Nz):
                                di = i1 - i
                                dj = j1 - j
                                dk = k1 - k
                                if di < 0:
                                    di += self.Nx
                                if dj < 0:
                                    dj += self.Ny
                                if dk < 0:
                                    dk += self.Nz
                                tmp[0] += 0.5 * self.Qxx[di, dj, dk] * dp[i1, j1, k1, 0]
                                tmp[1] += 0.5 * self.Qyy[di, dj, dk] * dp[i1, j1, k1, 1]
                                tmp[2] += 0.5 * self.Qzz[di, dj, dk] * dp[i1, j1, k1, 2]
                                tmp[1] += self.Qyz[di, dj, dk] * dp[i1, j1, k1, 2]
                                tmp[2] += self.Qzx[di, dj, dk] * dp[i1, j1, k1, 0]
                                tmp[0] += self.Qxy[di, dj, dk] * dp[i1, j1, k1, 1]
                    s += tmp @ dp[i, j, k, :]
        if normalize:
            s /= self.Nx * self.Ny * self.Nz
        return s

    def dp_int(self, dp: np.ndarray, normalize: bool = False) -> float:
        """Same as dp_int_direct but using FFT convolution (O(N log N))."""
        ux, uy, uz = dp[..., 0], dp[..., 1], dp[..., 2]
        ux_f = np.fft.fftn(ux)
        uy_f = np.fft.fftn(uy)
        uz_f = np.fft.fftn(uz)

        tmp0_f = (0.5 * self._Qxx_f.conj() * ux_f
                  + self._Qxy_f.conj() * uy_f)
        tmp1_f = (0.5 * self._Qyy_f.conj() * uy_f
                  + self._Qyz_f.conj() * uz_f)
        tmp2_f = (0.5 * self._Qzz_f.conj() * uz_f
                  + self._Qzx_f.conj() * ux_f)

        tmp0 = np.fft.ifftn(tmp0_f).real
        tmp1 = np.fft.ifftn(tmp1_f).real
        tmp2 = np.fft.ifftn(tmp2_f).real

        s = float(np.sum(tmp0 * ux + tmp1 * uy + tmp2 * uz))
        if normalize:
            s /= self.Nx * self.Ny * self.Nz
        return s


if __name__ == '__main__':
    import time

    qm = QMat('fesim/cache_q_10')
    dp = np.random.randn(qm.Nx, qm.Ny, qm.Nz, 3)

    t0 = time.perf_counter()
    slow = qm.dp_int_direct(dp)
    t1 = time.perf_counter()
    fast = qm.dp_int(dp)
    t2 = time.perf_counter()

    print(f'diff:        {abs(slow - fast):.2e}  (expected ~0)')
    print(f'direct O(N^2):   {t1 - t0:.4f} s')
    print(f'fft O(N log N):  {t2 - t1:.4f} s')
    print(f'speedup:         {(t1 - t0) / (t2 - t1):.1f}x')
