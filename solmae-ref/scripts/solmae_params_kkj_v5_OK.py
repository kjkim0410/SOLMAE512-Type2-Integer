# This Python file uses the following encoding: utf-8
"""
This script computes parameters and security estimates for SOLMAE.
 (Revised from Falcon's parameters.py. Kwangjo Kim. 23.8.2)
References:
- [ADPS16]: Post-quantum key exchange - a new hope
- [BDGL16]: ia.cr/2015/1128
- [DLP14]: ia.cr/2014/794
- [Duc18]:ia.cr/2017/999
- [EFG+22]: ia.cr/2021/1486, MITAKA
- [Falcon20]: https://falcon-sign.info
- [HPRR20]: ia.cr/2019/1411
- [Laa16]: https://pure.tue.nl/ws/files/14673128/20160216_Laarhoven.pdf
- [Lyu12]: ia.cr/2011/537
- [MR07]: https://cims.nyu.edu/~regev/papers/average.pdf
- [MW16]: ia.cr/2015/1123
- [NIST]: https://csrc.nist.gov/CSRC/media/Projects/Post-Quantum-Cryptography
          /documents/call-for-proposals-final-dec-2016.pdf
- [Pre17]: ia.cr/2017/480
- [SOL22]:https://solmae-sign.info
"""

from Crypto.Util.number import isPrime
from math import sqrt, exp, log, pi, floor, log
import sys

### HARDCODED CONSTANTS ########################
e = exp(1)
q = 1024 * 12 + 1 # modulo
NB_QUERIES = 2 ** 64 # NIST Recommendation
eta = 1/(2**41)
#############################################

DEBUG = True

# This is the maximal acceptable standard deviation for
# the individual Gaussians over Z, lifted from [HPRR20]
sigmax = 1.8205

# from parameters.py of Falcon python
def smooth(eps, n, normalized=True):
    """
    Compute the smoothing parameter eta_epsilon(Z^n).
    - if normalized is True, take the definition from [Pre17,Falcon]
    - if normalized is False, take the definition from [MR07]
    """
    rep = sqrt(log(2 * n * (1 + 1 / eps)) / pi)
    if normalized is True:
        return rep / sqrt(2 * pi)
    else:
        return rep

def ssmooth(eps, n):
    ''' Estimation of the somoothing parameter of ZZ^n '''
    return sqrt(log(2*n*(1+1/eps))/pi)/sqrt(2*pi)

##############################################################
# This comes from [Duc18] p.2
# by solving SVP in lattices of dim n while running 
# a sieve algorithm in projected sub-lattices of dim < n-d  
# where, by heuristic arguments   
# k \approx  (B * ln(4/3) / ln(B/(2*pi *exp))
# if B is 411, k= 37 this B'= B-k = 411-37 = 374, p.16 [Falcon20]
###################################################################

def dimensionsforfree(B):
    """
    d in [Duc18].
    """
    return round(B * log(4 / 3) / log(B / (2 * pi * exp(1))))

def print_security(B):
    '''
        Converts BIKZ security in bits for classical and quantum sieving
        using the core-svp methodology
    '''
    sec_qrec_classical  = floor(B*0.292)
    sec_qrec_quantum    = floor(B*0.265)
    print(" BIKZ:\t"+str(B))
    print(" Classical:\t"+str(sec_qrec_classical))
    print(" Quantum:\t"+str(sec_qrec_quantum))

input_para = [
{"d" :512, "alpha" : 1.17, "delta" : 0.065, "correction" : 1.04}
,
{"d" : 1024, "alpha" : 1.64, "delta" : 0.3,"correction": 1.04}
]

def compute_para(d,alp,delt,corr):
    gs_norm = alp * sqrt(q)
    smoothing = 1/pi * sqrt(1/2 * log(2*d * (1 + 1/eta)))
    sigma_sig = smoothing * alp * sqrt(q)
    gamma = corr * sigma_sig * sqrt(2*d)
    R_minus = (1/alp + delt) * sqrt(q)
    R_plus = (alp - delt) * sqrt(q)
    ## compute BKZ security using later by function 
    return gs_norm, smoothing, sigma_sig, gamma, R_minus, R_plus

for i in range(0, 2):
    d = input_para[i]["d"]   
    alp = input_para[i]["alpha"]
    delt = input_para[i]["delta"]
    corr = input_para[i]["correction"]

    gs_norm, smoothing, sigma_sig, gamma, R_minus, R_plus = compute_para(d,alp,delt,corr)

    if DEBUG:
        print("\n==== Security parameter for (d = {0}) \n".format(d))
        print("alpha(quality) = {: .3f}".format(alp))
        print("GS_norm = {: .3f}".format(gs_norm))
        print("smoothing = {: .3f}".format(smoothing))
        print("sigma_sig(sig width) = {: .3f}".format(sigma_sig))
        print("gamma = {: .3f}".format(gamma))
        print("gamma^2(sig. bound) = {: .3f}".format(gamma * gamma))
        print("R_minus = {: .3f}".format(R_minus))
        print("R_plus = {: .3f}".format(R_plus)) 

def solmae_security(n, sigma_offset, fg_norm, target_bitsec, \
        target_rejection = 0.1, verbose=True):
    '''
    SOLMAE parameter estimation
            - n is the dimension of the convolution ring (512 or 1024)
            - sigma_offset is the accessible standard deviation offset
            using the Hybrid sampler
            - fg_norm is the L_2 norm of the vector (f,g) in the canonical
              embedding
            - target_rejection is the maximal acceptable rate of rejection
            for signatures
    '''
    # Computation of the effective deviation of the hybrid sampler
    eps = 1/sqrt(target_bitsec * NB_QUERIES)
    sigma = sigma_offset*sqrt(q)*smooth(eps, 2*n)

    # Estimate of the signature size w.r.t rejection probability
    tau = 1.1
    while (1):
        max_sig_norm = floor(tau*sqrt(2*n)*sigma)
        rejection_rate = exp(2*n*(1-tau**2)/2)*tau**(2*n)
        if rejection_rate > target_rejection:   break
        else:                                   tau-= 0.001

    # Key recovery
    B = 100 # Initial blocksize
    sigma_fg = fg_norm/sqrt(2*n)
    while (1):
        left = (B/(2*pi*e))**(1-n/B) * sqrt(q)
        right = sqrt(3 * B/4)*sigma_fg
        if left > right: break
        else:            B += 1
    if verbose:
        print(" -----[ Key Recovery ] -----")
        print_security(B)

    # Signature forgery
    B = 100
    def condition_LH(beta):
        return min([(((pi*B)**(1/B)*B/(2 * pi * e)) ** ((2*n-k)/(2*B-2))) *
            q**(n/(2*n-k)) for k in range(n)])

    while condition_LH(B) > max_sig_norm: B += 1
    sec_forgery_classical = (B*0.292)
    sec_forgery_quantum   = (B*0.265)
    if verbose:
        print(" -----[ Signature forgery ] -----")
        print_security(B)
    return(sec_forgery_classical, sec_forgery_quantum)

print("\n == Security of SOLMAE-512==")
print((512,solmae_security(512, 2.04, 1.17*sqrt(q), 128, verbose=True)))

print("\n == Security of SOLMAE-1014==")
print((1024,solmae_security(1024, 2.33, 1.17*sqrt(q), 256, verbose=True)))