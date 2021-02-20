package integral;

import java.util.function.Function;

import utils.Constants;
import utils.SimpleMath;

/**
 * A numerical integration algorithm based on the Tanh-Sinh quadrature rule [1].
 * This implementation uses the method of weight an abcissa calculation in [2]
 * and further extensions in [3].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Takahasi, Hidetosi, and Masatake Mori. "Double exponential formulas
 * for numerical integration." Publications of the Research Institute for
 * Mathematical Sciences 9.3 (1974): 721-741.</li>
 * <li>[2] Mori, M.. “Developments in the Double Exponential Formulas for
 * Numerical Integration.” (1990).</li>
 * <li>[3] Bailey, D.. “Tanh-Sinh High-Precision Quadrature.” (2006).</li>
 * </ul>
 * </p>
 */
public final class TanhSinh extends Quadrature {

    private final int myLevels, myIters;
    private final double myRelTol;
    private final double[] myW, myX;

    /**
     * Creates a new instance of the Tanh-Sinh quadrature integrator. The weights
     * and abscissas are pre-computed and can be reused for evaluating multiple
     * integrals.
     * 
     * @param relativeTolerance the smallest acceptable relative change in integral
     *                          estimates in consecutive iterations that indicates
     *                          the algorithm has converged
     * @param absoluteTolerance the smallest acceptable absolute change in integral
     *                          estimates in consecutive iterations that indicates
     *                          the algorithm has converged
     * @param maxEvaluations    the maximum number of evaluations of each function
     *                          permitted
     * @param maxLevels         determines the number of interpolation points;
     *                          choosing large values may cause a heap exception;
     *                          defaults to 12
     */
    public TanhSinh(final double relativeTolerance, final double absoluteTolerance, final int maxEvaluations,
	    final int maxLevels) {
	super(absoluteTolerance, maxEvaluations);
	myRelTol = relativeTolerance;
	myLevels = maxLevels;

	// COMPUTE WEIGHTS AND ABSCISSAS
	final int len = ((1 << myLevels) << 2) * 5;
	myX = new double[len];
	myW = new double[len];
	myIters = tanhpts(2, Constants.EPSILON, myLevels, len, myW, myX);
    }

    /**
     * Creates a new instance of the Tanh-Sinh quadrature integrator. The weights
     * and abscissas are pre-computed and can be reused for evaluating multiple
     * integrals.
     * 
     * @param relativeTolerance the smallest acceptable relative change in integral
     *                          estimates in consecutive iterations that indicates
     *                          the algorithm has converged
     * @param absoluteTolerance the smallest acceptable absolute change in integral
     *                          estimates in consecutive iterations that indicates
     *                          the algorithm has converged
     * @param maxEvaluations    the maximum number of evaluations of each function
     *                          permitted
     */
    public TanhSinh(final double relativeTolerance, final double absoluteTolerance, final int maxEvaluations) {
	this(relativeTolerance, absoluteTolerance, maxEvaluations, 12);
    }

    public TanhSinh(final double tolerance, final int maxEvaluations) {
	this(100.0 * Constants.EPSILON, tolerance, maxEvaluations);
    }

    @Override
    final double properIntegral(final Function<? super Double, Double> f, final double a, final double b) {
	return tanhsinh(f, a, b, myW, myX);
    }

    @Override
    public final String getName() {
	return "TanhSinh";
    }

    private final double tanhsinh(final Function<? super Double, Double> f, final double a, final double b,
	    final double[] w, final double[] x) {

	// INITIALIZE CONSTANTS
	final double h = 0.5 * (b - a);
	final double mid = 0.5 * (a + b);

	// INITIALIZE ITERATION VARIABLES
	int fev = 0;
	double sum = 0.0;
	double estlast = 0.0;
	double v = 1.0;

	// MAIN LOOP OF TANH-SINH QUADRATURE
	for (int k = 1; k <= myLevels; ++k) {

	    // UPDATE ESTIMATE TO INCLUDE NEW LEVEL
	    final long stride = 1L << (myLevels - k);
	    for (int j = 0; j < myIters; j += stride) {
		if (k == 1 || (j & ((stride << 1) - 1)) != 0) {
		    if (j == 0) {
			final double fmid = f.apply(mid);
			sum += w[0] * fmid;
			++fev;
		    } else {
			double x1 = mid - h * x[j];
			double x2 = 2.0 * mid - x1;
			final double fx1 = f.apply(x1);
			final double fx2 = f.apply(x2);
			sum += w[j] * (fx1 + fx2);
			fev += 2;
		    }
		}
		if (fev >= myMaxEvals) {
		    myFEvals += fev;
		    return Double.NaN;
		}
	    }
	    v *= 0.5;
	    final double est = sum * h * v;

	    // CHECK CONVERGENCE
	    if (Double.isNaN(est)) {
		break;
	    }
	    final double tol = myRelTol * Math.abs(est) + myTol;
	    if (k > 1 && Math.abs(est - estlast) <= tol) {
		myFEvals += fev;
		return est;
	    }

	    // UPDATE LAST ESTIMATE
	    estlast = est;
	}

	// COULD NOT CONVERGE
	myFEvals += fev;
	return Double.NaN;
    }

    private static final int tanhpts(final int n, final double tol, final int levels, final int len, final double[] w,
	    final double[] x) {

	// INITIALIZE CONSTANTS
	final int hinv = 1 << levels;
	final double h = 1.0 / hinv;
	final double halfh = 0.5 * h;
	final double sinhh = Math.sinh(h);
	final double coshh = Math.cosh(h);
	final double sinhh2 = Math.sinh(halfh) * Math.PI;
	final double nest = Math.cosh(sinhh2 * Math.cosh(halfh));

	// INITIALIZE ITERATION VARIABLES
	double rlast = (coshh / nest) / nest;
	double kh = h;
	x[0] = 0.0;
	w[0] = Math.PI / 2.0;

	// MAIN LOOP OF WEIGHT AND ABSCISSA CALCULATION
	final double prec = SimpleMath.pow(tol, n);
	for (int k = 1; k < len; ++k) {

	    // RECURSION BY WATANABE
	    final double xi = Math.tanh(Math.PI / 2.0 * Math.sinh(kh));
	    final double s = sinhh2 * Math.cosh(kh + halfh);
	    final double den = Math.cosh(s) + xi * Math.sinh(s);
	    final double r = (coshh + sinhh * Math.tanh(kh)) / den / den;
	    x[k] = xi;
	    w[k] = w[k - 1] * rlast;

	    // TEST FOR CONVERGENCE
	    if (w[k] <= prec) {
		return k;
	    }

	    // UPDATE
	    rlast = r;
	    kh += h;
	}

	return len - 1;
    }
}
