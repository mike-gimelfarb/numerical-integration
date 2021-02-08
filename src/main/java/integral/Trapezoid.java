package integral;

import java.util.function.Function;

/**
 * A numerical integration method based on the trapezoid rule. Implements
 * iterative refinement to reduce the overall number of function evaluations.
 */
public final class Trapezoid extends Quadrature {

	public Trapezoid(final double tolerance) {
		super(tolerance);
	}

	@Override
	final double properIntegral(final Function<? super Double, Double> f, final double a, final double b) {
		return trapezoid_tol(f, a, b);
	}

	private double trapezoid_tol(final Function<? super Double, Double> f, final double a, final double b) {
		double q, qnew;
		int m = 1;
		q = 0.0;
		final int[] fev = new int[1];
		while (true) {
			qnew = trapezoid_refine(a, b, m, f, q, fev);
			if (1 < m && Math.abs(q - qnew) < myTol) {
				myFEvals += fev[0];
				return qnew;
			} else if (m >= 20) {
				myFEvals += fev[0];
				return Double.NaN;
			}
			q = qnew;
			++m;
		}
	}

	private static double trapezoid_refine(final double a, final double b, final int m,
			final Function<? super Double, Double> f, final double q, final int[] fev) {
		int i, i2, k, k2;
		double value, x;
		if (m < 1) {
			return Double.NaN;
		}
		if (m == 1) {
			value = (b - a) * 0.5 * (f.apply(a) + f.apply(b));
			fev[0] += 2;
		} else {
			k = 1 << (m - 2);
			k2 = k << 1;
			value = 0.0;
			for (i = 1; i <= k; ++i) {
				i2 = i << 1;
				x = ((k2 - i2 + 1) * a + (i2 - 1) * b) / k2;
				value += f.apply(x);
				++fev[0];
			}
			value = 0.5 * q + (b - a) * value / k2;
		}
		return value;
	}
}
