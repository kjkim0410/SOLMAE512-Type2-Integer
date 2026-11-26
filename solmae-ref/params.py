################################################################
# (Note) This is for SOLMAE purpose (23.8.10, kkj)
# SOLMAE_D is very special value to fix SOLMAE_512 or SOLMAE_1024
##################################################################

SOLMAE_D = 512
#SOLMAE_D = 1024

# Bytelength of the signing salt and header
HEAD_LEN = 1
SALT_LEN = 40  ## this must be 40 (Old :320)
SEED_LEN = 56

SOLMAE_Q = 12289

#Keep compatibility with Solmae C source
Params = {
    # SolmaeParam(512, 128)
    512: {
        "d": 512,
        "sigma": 1.32, #comes from Solmae C source
        "sigmin": 1.2778336969128337, #same as falcon need to check
        "sig_bound" : 33870790,
        "sig_bytelen": 666, #same as Falcon need to check
        "smoothing" : 1.338,
        "quality" : 1.17,
        "correction" : 0.065,
        "lower_radius": 101.95,
        "upper_radius" : 122.49,
        "signature_width" : 173.54,
        "slack" : 1.04,
    },
    # SolmaeParam(1024, 256)
    1024: {
        "d": 1024,
        "sigma": 1.32,
        "sigmin": 1.298280334344292, #same as Falcon need to check
        "sig_bound" : 134150669,
        "sig_bytelen": 1375, #same as Falcon need to check
        "smoothing" : 1.351,
        "quality" : 1.64,
        "correction" : 0.3,
        "lower_radius": 100.85,
        "upper_radius" : 148.54,
        "signature_width" : 245.62,
        "slack" : 1.04,
    },
}


logn = {
    2: 1,
    4: 2,
    8: 3,
    16: 4,
    32: 5,
    64: 6,
    128: 7,
    256: 8,
    512: 9,
    1024: 10
}
