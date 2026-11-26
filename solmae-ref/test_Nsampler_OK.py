from rng import ChaCha20
from os import urandom
from params import SOLMAE_D
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

m_ln2 = 0.693147180559945309417
m_pi = 3.14159265358979323846

def ffsll(word): # find first bit set in word and subtract 1. Program from math library of C source into python
    num = 0
    if word & 0xffffffff == 0:
        num += 32
        word >>= 32
    if word & 0xffff == 0:
        num += 16
        word >>= 16
    if word & 0xff ==0:
        num += 8
        word >>= 8
    if word & 0xf ==0:
        num += 4
        word >>= 4
    if word & 0x3 == 0:
        num += 2
        word >>= 2
    if word & 0x1 == 0:
        num += 1
    return num

def n_sampler():
    coeffs = []; seed_u = urandom(16)
    obj_u = ChaCha20(seed_u); seed_v = urandom(16)
    obj_v = ChaCha20(seed_v); seed_e = urandom(16)
    obj_e = ChaCha20(seed_e)

    for i in range(SOLMAE_D//2):
        # put 64bit random numbers in u, v, e.
        u = int.from_bytes(obj_u.randombytes(8), 'little')
        uf = 2* m_pi *(u & 0x1fffffffffffff) * 2**(-53)
        v = int.from_bytes(obj_v.randombytes(8), 'little')
        vf = 1/2+ (v & 0x1fffffffffffff)* 2**(-54)
        e_1 = int.from_bytes(obj_e.randombytes(8), 'little')
        e_2 = int.from_bytes(obj_e.randombytes(8), 'little')
        # Use constant time macro
        x = 64 + ffsll(e_2); y = ffsll(e_1)
        c = ((~e_1)&(e_1 - 1)) >> 63
        geom = (x & (-(c&1)))^(y & (~(-(c&1))))
        vf = np.sqrt(SOLMAE_D*(m_ln2*geom-np.log(vf)))
        # statement above is equivalent to u_rho ~ U(0,1), vf = sqrt(-d*ln(u_rho))
        coeffs.append(vf*np.cos(uf)); coeffs.append(vf*np.sin(uf))
    return coeffs

if __name__ == '__main__':
    # generate datas
    x_list = []; y_list = []; datas = []
    n = 1000 # Set the number of test 
    for _ in range(n):
        coeffs = n_sampler(); data = []
        for i in range(SOLMAE_D//2):
            x = coeffs[2*i]; y = coeffs[2*i+1]
            x_list.append(x); y_list.append(y)
            data.append(x);data.append(y)
        datas.append(data)

    # plot the data in 2D
    fig = plt.figure(figsize = (20, 40)); ax1 = fig.add_subplot(2, 1, 1)
    plt.title("Distribution represented in 2D")
    plt.plot(x_list,y_list, 'o', markersize=3)
    plt.xlim([-80, 80]); plt.ylim([-80, 80])

    # calculate sample mean and sample variance
    datas = np.array(datas);  mu = np.array([0 for _ in range(SOLMAE_D)])
    cov_matrix = SOLMAE_D/2*np.eye(SOLMAE_D)
    cov_matrix_inv = np.linalg.inv(cov_matrix)
    
    data_mean = np.mean(datas, axis = 0)
    print("mean:", data_mean)

    data_cov_matrix = np.cov(datas.T, ddof = 1)
    print("covariance: ", data_cov_matrix)
    data_cov_matrix_inv = np.linalg.inv(data_cov_matrix)

    # check if the datas follow the distribution N_{d/2} with qq plot
    ax2 = fig.add_subplot(2, 1, 2)
    Maha_dist = []
    
    for data in datas:
        data = np.array(data)
        Maha_dist.append(float(np.dot(np.dot(data.T, cov_matrix_inv), data)))
    Maha_dist = np.array(Maha_dist)
    stats.probplot(Maha_dist, dist = stats.chi2(SOLMAE_D), plot=ax2)
    plt.title("chi-square QQ-plot")
    plt.show()