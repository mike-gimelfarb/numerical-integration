package integral;

import java.util.function.Function;

import utils.Constants;

/**
 * Implements an adaptive numerical integrator based on the 7-point
 * Gauss-Lobatto rule, as introduced in [1].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Gander, W., Gautschi, W. Adaptive Quadrature—Revisited. BIT Numerical
 * Mathematics 40, 84–101 (2000). https://doi.org/10.1023/A:1022318402393</li>
 * </ul>
 * </p>
 */
public final class GaussLobatto extends Quadrature {

    private static final double ALPHA = Constants.SQRT2 / Constants.SQRT3;
    private static final double BETA = 1.0 / Constants.SQRT5;

    private static final double[] X = { 0.94288241569547971906, 0.64185334234578130578, 0.23638319966214988028 };
    private static final double[] Y = { 0.0158271919734801831, 0.0942738402188500455, 0.1550719873365853963,
	    0.1888215739601824544, 0.1997734052268585268, 0.2249264653333395270, 0.2426110719014077338 };
    private static final double[] C = { 77.0, 432.0, 625.0, 672.0 };

    public GaussLobatto(final double tolerance) {
	super(tolerance);
    }

    @Override
    final double properIntegral(final Function<? super Double, Double> f, final double a, final double b) {
	return dlob8e(f, a, b);
    }

    private double dlob8e(final Function<? super Double, Double> f, final double a, final double b) {

	// COMPUTE INTERPOLATION POINTS
	final double mid = 0.5 * (a + b);
	final double h = 0.5 * (b - a);
	final double y1 = f.apply(a);
	final double y3 = f.apply(mid - h * ALPHA);
	final double y5 = f.apply(mid - h * BETA);
	final double y7 = f.apply(mid);
	final double y9 = f.apply(mid + h * BETA);
	final double y11 = f.apply(mid + h * ALPHA);
	final double y13 = f.apply(b);
	final double f1 = f.apply(mid - h * X[0]);
	final double f2 = f.apply(mid + h * X[0]);
	final double f3 = f.apply(mid - h * X[1]);
	final double f4 = f.apply(mid + h * X[1]);
	final double f5 = f.apply(mid - h * X[2]);
	final double f6 = f.apply(mid + h * X[2]);
	final double est1 = (y1 + y13 + 5.0 * (y5 + y9)) * (h / 6.0);
	final double est2 = (C[0] * (y1 + y13) + C[1] * (y3 + y11) + C[2] * (y5 + y9) + C[3] * y7) * (h / 1470.0);
	myFEvals += 13;

	// COMPUTE ERROR ESTIMATE
	final double s = (Y[0] * (y1 + y13) + Y[1] * (f1 + f2) + Y[2] * (y3 + y11) + Y[3] * (f3 + f4) + Y[4] * (y5 + y9)
		+ Y[5] * (f5 + f6) + Y[6] * y7) * h;
	final double r = Math.abs(est1 - s) != 0.0 ? Math.abs(est2 - s) / Math.abs(est1 - s) : 1.0;
	final double rtol = r > 0.0 && r < 1.0 ? myTol / r : myTol;
	final double s1 = s == 0.0 ? Math.abs(b - a) : s;

	// CALL MAIN SUBROUTINE
	return dlob8(f, a, b, y1, y13, s1, rtol);
    }

    private double dlob8(final Function<? super Double, Double> f, final double a, final double b, final double fa,
	    final double fb, final double s, final double rtol) {

	// INITIALIZE
	final double h = 0.5 * (b - a);
	final double mid = 0.5 * (a + b);
	final double mll = mid - h * ALPHA;
	final double ml = mid - h * BETA;
	final double mr = mid + BETA * h;
	final double mrr = mid + h * ALPHA;
	final double fmll = f.apply(mll);
	final double fml = f.apply(ml);
	final double fmid = f.apply(mid);
	final double fmr = f.apply(mr);
	final double fmrr = f.apply(mrr);
	final double est2 = (C[0] * (fa + fb) + C[1] * (fmll + fmrr) + C[2] * (fml + fmr) + C[3] * fmid) * (h / 1470.0);
	final double est1 = (fa + fb + 5.0 * (fml + fmr)) * (h / 6.0);
	myFEvals += 8;

	// CHECK CONVERGENCE
	if (Math.abs(est1 - est2) <= myTol * s || mll <= a || b <= mrr) {
	    return est2;
	}

	// RECURSION
	return dlob8(f, a, mll, fa, fmll, s, rtol) + dlob8(f, mll, ml, fmll, fml, s, rtol)
		+ dlob8(f, ml, mid, fml, fmid, s, rtol) + dlob8(f, mid, mr, fmid, fmr, s, rtol)
		+ dlob8(f, mr, mrr, fmr, fmrr, s, rtol) + dlob8(f, mrr, b, fmrr, fb, s, rtol);
    }
}
