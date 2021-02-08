package series.dp;

import utils.SimpleMath;

/**
 * Implements the family of Generalized Levin algorithms for convergence of
 * infinite series based on remainder sequences described in [1, 2] and as
 * outlined in [3].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Levin, David. "Development of non-linear transformations for
 * improving convergence of sequences." International Journal of Computer
 * Mathematics 3.1-4 (1972): 371-388.</li>
 * <li>[2] Smith, David A., and William F. Ford. "Acceleration of linear and
 * logarithmic convergence." SIAM Journal on Numerical Analysis 16.2 (1979):
 * 223-240.</li>
 * <li>[3] Weniger, Ernst Joachim. "Nonlinear sequence transformations for the
 * acceleration of convergence and the summation of divergent series." arXiv
 * preprint math/0306302 (2003).</li>
 * </ul>
 * </p>
 */
public final class Levin extends SeriesAlgorithm {

	/**
	 * The type of remainder sequence to assume when accelerating the convergence of
	 * a series.
	 */
	public static enum RemainderSequence {

		T {
			@Override
			public double remainder(int n, double a1, double a2, double beta) {
				return a2;
			}
		},
		U {
			@Override
			public double remainder(int n, double a1, double a2, double beta) {
				return a1 * (beta + n);
			}
		},
		V {
			@Override
			public double remainder(int n, double a1, double a2, double beta) {
				return a1 * a2 / (a1 - a2);
			}
		};

		public abstract double remainder(int n, double a1, double a2, double beta);
	}

	private final String prefix;
	private final double myBeta;
	private final RemainderSequence myRemainder;

	private double myE, myTerm;
	private final double[] myTabHi, myTabLo;

	/**
	 * Creates a new instance of the Generalized Levin algorithm with customizable
	 * remainder sequences.
	 * 
	 * @param tolerance the smallest acceptable change in series evaluation in
	 *                  consecutive iterations that indicates the algorithm has
	 *                  converged
	 * @param maxIters  the maximum number of sequence terms to evaluate before
	 *                  giving up
	 * @param beta      the constant described in [1, 3] for the Levin u-transform
	 *                  remainder sequence: if this sequence is not selected this
	 *                  argument is ignored
	 * @param remainder the type of assumed remainder sequence as described in [3]
	 */
	public Levin(final double tolerance, final int maxIters, final double beta, final RemainderSequence remainder) {
		super(tolerance, maxIters);
		myBeta = beta;
		myRemainder = remainder;
		prefix = remainder.name();

		myTabHi = new double[maxIters];
		myTabLo = new double[maxIters];
	}

	/**
	 * Creates a new instance of the Generalized Levin algorithm with customizable
	 * remainder sequences.
	 * 
	 * @param tolerance the smallest acceptable change in series evaluation in
	 *                  consecutive iterations that indicates the algorithm has
	 *                  converged
	 * @param maxIters  the maximum number of sequence terms to evaluate before
	 *                  giving up
	 * @param remainder the type of assumed remainder sequence as described in [3]
	 */
	public Levin(final double tolerance, final int maxIters, final RemainderSequence remainder) {
		this(tolerance, maxIters, 1.0, remainder);
	}

	@Override
	public final double next(final double e, final double term) {

		// shift over elements by one by skipping first iteration so that
		// the iteration is actually computing the term one over from the
		// one used in the current computation
		if (myIndex == 0) {
			myE = e;
			myTerm = term;
			++myIndex;
			return Double.NaN;
		}

		// compute the remainder term
		final int k = myIndex - 1;
		final double rem = myRemainder.remainder(k + 1, myE, e, myBeta);

		// update table
		myTabHi[k] = myTerm / rem;
		myTabLo[k] = 1.0 / rem;
		if (k > 0) {
			final int itm1 = k - 1;
			myTabHi[itm1] = myTabHi[k] - myTabHi[itm1];
			myTabLo[itm1] = myTabLo[k] - myTabLo[itm1];
			if (k > 1) {
				final double bn1 = myBeta + k - 1.0;
				final double bn2 = bn1 + 1.0;
				final double coeff = bn1 / bn2;
				for (int j = 2; j <= k; ++j) {
					final int itmj = k - j;
					final double cjm2 = SimpleMath.pow(coeff, j - 2);
					final double fac = (myBeta + itmj) * cjm2 / bn2;
					myTabHi[itmj] = myTabHi[itmj + 1] - fac * myTabHi[itmj];
					myTabLo[itmj] = myTabLo[itmj + 1] - fac * myTabLo[itmj];
				}
			}
		}

		// increment counters and temporaries
		myE = e;
		myTerm = term;
		++myIndex;

		// correct for underflow
		if (Math.abs(myTabLo[0]) < TINY) {
			return HUGE;
		} else {
			return myTabHi[0] / myTabLo[0];
		}
	}

	@Override
	public String getName() {
		return "Levin" + " " + prefix;
	}
}
