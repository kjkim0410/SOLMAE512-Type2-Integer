from fxp import Fxp, add, sub, mul, div
from fft_fxp import fft_fxp, ifft_fxp
from fxp_formats import F


def main():
    a = Fxp.from_float(1.25, 30)
    b = Fxp.from_float(-0.75, 30)
    assert abs(add(a, b, 30).to_float() - 0.5) < 1e-8
    assert abs(sub(a, b, 30).to_float() - 2.0) < 1e-8
    assert abs(mul(a, b, 30).to_float() + 0.9375) < 1e-8
    assert abs(div(a, Fxp.from_float(0.5, 30), 30).to_float() - 2.5) < 1e-8

    x = [1, -2, 3, 4, -1, 0, 2, -3]
    y = [z.to_float() for z in ifft_fxp(fft_fxp(x, F["TARGET"]), F["COEFF"])]
    assert max(abs(u-v) for u, v in zip(x, y)) < 1e-4
    print("Type-2 FxP core test: PASS")


if __name__ == "__main__":
    main()
