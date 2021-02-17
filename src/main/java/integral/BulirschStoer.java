package integral;

import java.util.Arrays;
import java.util.Stack;
import java.util.function.Function;

import utils.Constants;

/**
 * An adaptive integration routine that estimates the integral on each region
 * using the trapezoid rule for different step sizes, then a rational function
 * is fit to the result. The integral on the subinterval is estimated by
 * extrapolating when the step size goes to zero [1]. The subroutine iratex is
 * translated from Fortran code by John Burkardt.
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Bulirsch, R., Stoer, J.: Fehlerabschätzungen und Extrapolation mit
 * rationalen Funktionen bei Verfahren vom Richardson-Typus. Numer. Math.6,
 * 413---427 (1964)</li>
 * </ul>
 * </p>
 */
public final class BulirschStoer extends Quadrature {

    private final int myMaxEvals;

    /**
     * Creates a new instance of the Bulirsch-Stoer quadrature integrator.
     * 
     * @param tolerance      the smallest acceptable absolute change in integral
     *                       estimates in consecutive iterations that indicates the
     *                       algorithm has converged
     * @param maxEvaluations maximum number of functions evaluations
     */
    public BulirschStoer(final double tolerance, final int maxEvaluations) {
	super(tolerance);
	myMaxEvals = maxEvaluations;
    }

    public BulirschStoer(final double tolerance) {
	this(tolerance, 9999);
    }

    @Override
    final double properIntegral(final Function<? super Double, Double> f, final double a, final double b) {
	return bulirschStoer(f, a, b);
    }

    private double bulirschStoer(final Function<? super Double, Double> func, final double a, final double b) {

	// keep a stack of interval partitions which need to be divided
	final Stack<double[]> stack = new Stack<>();
	stack.push(new double[] { a, b });

	double result = 0.0;
	final int[] fev = new int[1];

	// keep arrays in which to store estimate of the integral and error
	final double[] estout = new double[1];
	final double[] epsout = new double[1];
	final int[] indout = new int[1];

	// work arrays
	final double[] work1 = new double[6];
	final double[] work2 = new double[7];

	// we adapt the algorithm by Bulirsch and Stoer which provides
	// error estimates into a globally adaptive method
	while (!stack.isEmpty()) {

	    // remove the top interval from stack
	    final double[] interval = stack.pop();
	    final double a1 = interval[0];
	    final double b1 = interval[1];

	    // estimate the integral and error over this interval
	    final double epsin = (b1 - a1) * myTol;
	    iratex(func, a1, b1, epsin, epsout, estout, indout, fev, work1, work2);

	    // empty intervals are ignored
	    if (b1 - a1 == 0.0) {
		continue;
	    }

	    // check if the accuracy was reached
	    if (Math.abs(epsout[0]) <= epsin && indout[0] == 1) {
		result += estout[0];
	    } else {

		// subdivide the interval
		final double mid = 0.5 * (a1 + b1);
		stack.push(new double[] { a1, mid });
		stack.push(new double[] { mid, b1 });
	    }

	    // check budget
	    if (fev[0] >= myMaxEvals) {
		myFEvals += fev[0];
		return Double.NaN;
	    }
	}
	myFEvals += fev[0];
	return result;
    }

    private static void iratex(final Function<? super Double, Double> func, final double a, final double b,
	    double epsin, final double[] epsout, final double[] result, final int[] ind, final int[] fev,
	    final double[] d, final double[] dt) {
	double arg, ba, c, dl, ddt, den, e, ent, gr, hm, rnderr, sm, t, t1, t2, t2a, ta, tab = 0.0, tb, tnt, v = 0.0, w;
	int i, m, mr, n, np1;
	boolean bo, bu, odd;

	if (a == b) {
	    result[0] = 0.0;
	    return;
	}

	Arrays.fill(d, 0.0);
	Arrays.fill(dt, 0.0);
	rnderr = Constants.EPSILON;
	epsin = Math.max(epsin, 8.0 * rnderr);
	ind[0] = 0;
	n = 2;
	np1 = 3;
	ba = b - a;
	t1 = gr = sm = 0.0;
	t2a = 0.5 * (func.apply(a) + func.apply(b));
	fev[0] += 2;
	t2 = t2a;
	tb = Math.abs(t2a);
	c = t2 * ba;
	dt[1 - 1] = c;
	odd = true;
	bu = false;

	for (m = 1; m <= 15; ++m) {
	    bo = 7 <= m;
	    hm = ba / n;

	    // N+1 is odd.
	    if (odd) {
		for (i = 1; i <= n; i += 2) {
		    arg = a + i * hm;
		    w = func.apply(arg);
		    ++fev[0];
		    t2 += w;
		    tb += Math.abs(w);
		}
		ent = t2;
		tab = tb * Math.abs(hm);
		d[1 - 1] = 16.0 / 9.0;
		d[3 - 1] = 4.0 * d[1 - 1];
		d[5 - 1] = 4.0 * d[3 - 1];
	    } else {

		// N+1 is even.
		for (i = 1; i <= n; i += 6) {
		    w = i * hm;
		    t1 += func.apply(a + w) + func.apply(b - w);
		    fev[0] += 2;
		}
		ent = t1 + t2a;
		t2a = t2;
		d[1 - 1] = 2.25;
		d[3 - 1] = 9.0;
		d[5 - 1] = 36.0;
	    }

	    ddt = dt[1 - 1];
	    t = ent * hm;
	    dt[1 - 1] = ent = t;
	    if (m < 7) {
		mr = m;
		w = n * n;
		d[m - 1] = w;
	    } else {
		mr = 6;
		d[6 - 1] = 64.0;
		w = 144.0;
	    }

	    for (i = 1; i <= mr; ++i) {
		dl = d[i - 1] * ddt;
		den = dl - ent;
		e = ent - ddt;
		tnt = ent;
		v = ent = 0.0;
		if (epsin <= Math.abs(den)) {
		    e /= den;
		    v = tnt * e;
		    ent = dl * e;
		    t += v;
		    ddt = dt[i + 1 - 1];
		}
		dt[i + 1 - 1] = v;
	    }

	    ta = c;
	    c = t;
	    result[0] = c;
	    if (!bo) {
		t -= v;
	    }
	    v = t - ta;
	    t += v;
	    epsout[0] = Math.abs(v);
	    if (ta < t) {
		dl = ta;
		ta = t;
		t = dl;
	    }
	    bo = bo || (ta < gr && sm < t);
	    if (bu && bo && epsout[0] < epsin * w * tab) {
		v = rnderr * tab;
		epsout[0] = Math.max(epsout[0], v);
		if (bo) {
		    ind[0] = 1;
		}
		return;
	    }
	    gr = ta + epsin;
	    sm = t - epsin;
	    odd = !odd;
	    n = np1;
	    np1 = n + 1;
	    bu = bo;
	    d[2 - 1] = 4.0;
	    d[4 - 1] = 16.0;
	}
	bo = false;
	v = rnderr * tab;
	epsout[0] = Math.max(epsout[0], v);
	if (bo) {
	    ind[0] = 1;
	}
    }
}
