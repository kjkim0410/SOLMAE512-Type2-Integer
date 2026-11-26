from math import sqrt, exp, log, pi, floor

# Constants
e = exp(1)
q = 1024 * 12 + 1  # Modulo
NB_QUERIES = 2 ** 64  # NIST Recommendation

eta = 1 / (2 ** 41)
sigmax = 1.8205  # Max. acceptable std. dev for Gaussian sampler

DEBUG = True

def smooth(eps, n, normalized=True):
    """Calculate the smoothing parameter."""
    rep = sqrt(log(2 * n * (1 + 1 / eps)) / pi)
    if normalized:
        return rep / sqrt(2 * pi)
    else:
        return rep

def ssmooth(eps, n):
    """Estimation of the smoothing parameter of ZZ^n."""
    return sqrt(log(2 * n * (1 + 1 / eps)) / pi) / sqrt(2 * pi)

def dimensionsforfree(B):
    """Calculate dimensions for free."""
    return round(B * log(4 / 3) / log(B / (2 * pi * exp(1))))

def print_security(B):
    """Print the security parameters."""
    sec_qrec_classical = floor(B * 0.292)
    sec_qrec_quantum = floor(B * 0.265)
    print(f"BIKZ:\t{B}")
    print(f"Classical:\t{sec_qrec_classical}")
    print(f"Quantum:\t{sec_qrec_quantum}")

def compute_para(d, alp, delt, corr):
    """Compute parameters based on input."""
    gs_norm = alp * sqrt(q)
    smoothing = 1 / pi * sqrt(1 / 2 * log(2 * d * (1 + 1 / eta)))
    sigma_sig = smoothing * alp * sqrt(q)
    gamma = corr * sigma_sig * sqrt(2 * d)
    R_minus = (1 / alp + delt) * sqrt(q)
    R_plus = (alp - delt) * sqrt(q)
    return gs_norm, smoothing, sigma_sig, gamma, R_minus, R_plus

def solmae_security(n, sigma_offset, fg_norm, target_bitsec, target_rejection=0.1, verbose=True):
    """Calculate security parameters for SOLMAE."""
    eps = 1 / sqrt(target_bitsec * NB_QUERIES)
    sigma = sigma_offset * sqrt(q) * smooth(eps, 2 * n)
    tau = 1.1  # Estimate of the signature size w.r.t rejection prob.
    
    while True:
        max_sig_norm = floor(tau * sqrt(2 * n) * sigma)
        rejection_rate = exp(2 * n * (1 - tau ** 2) / 2) * tau ** (2 * n)
        if rejection_rate > target_rejection:
            break
        else:
            tau -= 0.001
    
    B = 100  # Initial block for Key recovery
    sigma_fg = fg_norm / sqrt(2 * n)
    
    while True:
        left = (B / (2 * pi * e)) ** (1 - n / B) * sqrt(q)
        right = sqrt(3 * B / 4) * sigma_fg
        if left > right:
            break
        else:
            B += 1
    
    if verbose:
        print(" -----[ Key Recovery ] -----")
        print_security(B)
    
    B = 100  # Signature forgery
    
    def condition_LH(beta):
        return min([(((pi * B) ** (1 / B) * B / (2 * pi * e)) ** ((2 * n - k) / (2 * B - 2))) *
                    q ** (n / (2 * n - k)) for k in range(n)])
    
    while condition_LH(B) > max_sig_norm:
        B += 1
    
    sec_forgery_classical = (B * 0.292)
    sec_forgery_quantum = (B * 0.265)
    
    if verbose:
        print(" -----[ Signature forgery ] -----")
        print_security(B)
    return 

def para_gen():
    input_para = [
        {"d": 512, "alpha": 1.17, "delta": 0.065, "correction": 1.04},
        {"d": 1024, "alpha": 1.64, "delta": 0.3, "correction": 1.04}
    ]
    
    for params in input_para:
        d = params["d"]
        alp = params["alpha"]
        delt = params["delta"]
        corr = params["correction"]
        gs_norm, smoothing, sigma_sig, gamma, R_minus, R_plus = compute_para(d, alp, delt, corr)
    
        if DEBUG:
            print(f"\n ** Security parameters for d = {d} **")
            print("alpha(quality)     = {:.3f}".format(alp).ljust(25))
            print("GS_norm            = {:.3f}".format(gs_norm).ljust(25))
            print("smoothing          = {:.3f}".format(smoothing).ljust(25))
            print("sigma_sig(sig width) = {:.3f}".format(sigma_sig).ljust(25))
            print("gamma              = {:.3f}".format(gamma).ljust(25))
            print("gamma^2(sig. bound) = {:.3f}".format(gamma * gamma).ljust(25))
            print("R_minus            = {:.3f}".format(R_minus).ljust(25))
            print("R_plus             = {:.3f}".format(R_plus).ljust(25))
       
    print("\n== C/Q security of SOLMAE_512 ==")
    solmae_security(512, 2.04, 1.17 * sqrt(q), 128, verbose=True)
    
    print("\n== C/Q security of SOLMAE_1014 ==")
    solmae_security(1024, 2.33, 1.17 * sqrt(q), 256, verbose=True)

if __name__ == "__main__":
    para_gen()

print("Test passed.")