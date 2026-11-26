#############################################################
#This is to generate all keys and its test program for SOLMAE.
#############################################################
from Pairgen import Pairgen
from ntt import ntt, intt, div_ntt
from ntrugen import ntru_solve
from fft import fft, ifft, add_fft, sub_fft, mul_fft, div_fft, adj_fft, cut_half_fft
from numpy import sqrt
from params import Params, SOLMAE_D
from os import urandom

class secret_key:
    def __init__(self, n=None, fgFG=None):
        self.f = []; self.g = []; self.F = []; self.G = []
        self.Sigma1 = []; self.Sigma2 = []
        self.b10_fft = []; self.b11_fft = []; self.b20_fft = []; self.b21_fft = []
        self.beta10_fft = []; self.beta11_fft = []
        self.beta20_fft = []; self.beta21_fft = []
  
        if fgFG is not None:
            self.f, self.g, self.F, self.G = fgFG

    # def __repr__(self):
    #     res = "Secretkey\n======================\n"
    #     res += f"f      = {self.f}\n======================\n"
    #     res += f"g      = {self.g}\n======================\n"
    #     res += f"F      = {self.F}\n======================\n"
    #     res += f"G      = {self.G}\n======================\n"
    #     res += f"Sigma1 = {self.Sigma1}\n======================\n"
    #     res += f"Sigma2 = {self.Sigma2}\n======================\n"
    #     res += f"b10_fft = {self.b10_fft}\n======================\n"
    #     res += f"b11_fft = {self.b11_fft}\n======================\n"
    #     res += f"b20_fft = {self.b20_fft}\n======================\n"
    #     res += f"b21_fft = {self.b21_fft}\n======================\n"
    #     res += f"beta10_fft = {self.beta10_fft}\n=======================\n"
    #     res += f"beta11_fft = {self.beta11_fft}\n======================\n"
    #     res += f"beta20_fft = {self.beta20_fft}\n======================\n"
    #     res += f"beta21_fft = {self.beta21_fft}\n======================\n"
    #     return res

class public_key:
    def __init__(self):
        self.h = []
    
    # def __repr__(self):
    #     return f"===============\npublickey h = {self.h}\n=================\n"

def keygen(randombytes=urandom):
    # sk = secret_key(f,g,F,G)  # Error
    sk = secret_key()
    pk = public_key()
    while True:
        f, g, _ = Pairgen(randombytes)
        try:
            f_ntt = ntt(f)
            g_ntt = ntt(g)
            h_ntt = div_ntt(g_ntt, f_ntt)
        except ZeroDivisionError:
            continue
        try:
            F, G = ntru_solve(f, g)
        except ValueError:
            continue

        sk.f = f; sk.g = g; sk.F = F; sk.G = G  
        pk.h = intt(h_ntt)
        break

    # Consistency values (2023.8.10, kkj)
    eta_sq = Params[SOLMAE_D]["smoothing"] ** 2
    sig_width = Params[SOLMAE_D]["signature_width"] ** 2
    eta_sq_fft = [eta_sq for _ in range(SOLMAE_D)]
    sig_width_fft = [sig_width for _ in range(SOLMAE_D)]
    
    sk.b10_fft = fft(f); sk.b11_fft = fft(g)
    
    b1_norm = add_fft(mul_fft(adj_fft(sk.b10_fft), sk.b10_fft), mul_fft(adj_fft(sk.b11_fft), sk.b11_fft))
    sk.beta10_fft = div_fft(sk.b10_fft, b1_norm)
    sk.beta11_fft = div_fft(sk.b11_fft, b1_norm)
    sk.Sigma1 = [sqrt(elem) for elem in cut_half_fft(sub_fft(div_fft(sig_width_fft, b1_norm), eta_sq_fft))]

    sk.b20_fft = fft(F); sk.b21_fft = fft(G)
    
    temp_fft = add_fft(mul_fft(adj_fft(sk.beta10_fft), sk.b20_fft), mul_fft(adj_fft(sk.beta11_fft), sk.b21_fft))
    sk.b20_tild_fft = sub_fft(sk.b20_fft, mul_fft(temp_fft, sk.b10_fft))
    sk.b21_tild_fft = sub_fft(sk.b21_fft, mul_fft(temp_fft, sk.b11_fft))
    
    b2_tild_norm = add_fft(mul_fft(adj_fft(sk.b20_tild_fft), sk.b20_tild_fft), mul_fft(adj_fft(sk.b21_tild_fft), sk.b21_tild_fft))
    sk.beta21_fft = div_fft(sk.b21_tild_fft, b2_tild_norm)
    sk.Sigma2 = [sqrt(elem) for elem in cut_half_fft(sub_fft(div_fft(sig_width_fft, b2_tild_norm), eta_sq_fft))]

    return sk, pk

if __name__ == '__main__':
    sk, pk = keygen()
    print("== Leading 3 to 10 values of key pairs for SOLMAE_",SOLMAE_D,"==")
    print("f            = ".ljust(8) + "[" + ", ".join(["{}".format(x) for x in sk.f[:10]]) + ", ...]") 
    print("g            = ".ljust(8) + "[" + ", ".join(["{}".format(x) for x in sk.g[:10]]) + ", ...]")
    print("F            = ".ljust(8) + "[" + ", ".join(["{}".format(x) for x in sk.F[:10]]) + ", ...]") 
    print("G            = ".ljust(8) + "[" + ", ".join(["{}".format(x) for x in sk.G[:10]]) + ", ...]")   
    print("h            = ".ljust(8) + "[" + ", ".join(["{}".format(x) for x in pk.h[:9]]) + ", ...]") 
    print("f_fft        = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b10_fft[:3]]) + ", ...]")
    print("g_fft        = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b11_fft[:3]]) + ", ...]")
    print("F_fft        = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b20_fft[:3]]) + ", ...]")
    print("G_fft        = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b21_fft[:3]]) + ", ...]")
    print("beta10_fft   = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b10_fft[:3]]) + ",  ...]")
    print("beta11_fft   = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b11_fft[:3]]) + ", ...]")
    print("beta20_fft   = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b20_fft[:3]]) + ", ...]")
    print("beta21_fft   = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.b21_fft[:3]]) + ", ...]")
    print("Sigma1       = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.Sigma1[:3]]) + ", ...]")
    print("Sigma2       = ".ljust(8) + "[" + ", ".join(["{:.3f}".format(x) for x in sk.Sigma2[:3]]) + ", ...]")
