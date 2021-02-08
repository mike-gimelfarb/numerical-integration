package series.dp;

/**
 * Implements the Wynn Epsilon algorithm for convergence of infinite series
 * introduced in [1] and using an efficient moving lozenge technique as
 * described in [2].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Wynn, Peter. "On a device for computing the e m (S n)
 * transformation." Mathematical Tables and Other Aids to Computation (1956):
 * 91-96.</li>
 * <li>[2] Weniger, Ernst Joachim. "Nonlinear sequence transformations for the
 * acceleration of convergence and the summation of divergent series." arXiv
 * preprint math/0306302 (2003).</li>
 * </ul>
 * </p>
 */
public final class WynnEpsilon extends SeriesAlgorithm {

	private final double[] myTab;

	public WynnEpsilon(final double tolerance, final int maxIters) {
		super(tolerance, maxIters);
		myTab = new double[maxIters];
	}

	@Override
	public final double next(final double e, final double term) {

		// initialization
		myTab[myIndex] = term;
		if (myIndex == 0) {
			++myIndex;
			return term;
		}

		// main loop of Wynn Epsilon algorithm
		double aux = 0.0;
		for (int j = myIndex; j >= 1; --j) {

			// compute the next epsilon
			final double temp = aux;
			aux = myTab[j - 1];
			double diff = myTab[j] - aux;

			// correct denominators for underflow
			if (Math.abs(diff) <= TINY) {
				diff = HUGE;
			} else {
				diff = temp + 1.0 / diff;
			}
			myTab[j - 1] = diff;
		}

		// prepare result
		final double result = myTab[myIndex & 1];
		++myIndex;
		return result;
	}

	@Override
	public String getName() {
		return "Wynn Epsilon";
	}
}
