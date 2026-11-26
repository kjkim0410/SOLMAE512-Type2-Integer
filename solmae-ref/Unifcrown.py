###########################
# This is for Solmae only.
###########################
import numpy as np
from rng import ChaCha20
from os import urandom
import matplotlib.pyplot as plt

def Unifcrown(R_min, R_max, randombytes=urandom):
    u_rho = int.from_bytes(randombytes(8), 'little')
    u_theta = int.from_bytes(randombytes(8), 'little')
    u_rho = (u_rho & 0x1fffffffffffff) * 2**(-53)
    u_theta = (u_theta & 0x1fffffffffffff) * 2**(-53)
    
    rho = np.sqrt(R_min**2+u_rho*(R_max**2-R_min**2))
    x = rho*np.cos(np.pi/2*u_theta)
    y = rho*np.sin(np.pi/2*u_theta) # equivalent to Algorithm 9
    return x, y

if __name__ == '__main__':
    x_list = []
    y_list = []
    for _ in range(5000):
        x, y = Unifcrown(2, 5)
        x_list.append(x)
        y_list.append(y)

    plt.plot(x_list,y_list, 'o', markersize=3)

    plt.show()