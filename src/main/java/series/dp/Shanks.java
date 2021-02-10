package series.dp;

/**
 * Implements the Shanks transform for accelerating convergence of infinite
 * series. Implemented using the G-transform as described in [1].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Varga, R.S. Extrapolation methods: theory and practice. Numer Algor
 * 4, 305 (1993). https://doi.org/10.1007/BF02144109</li>
 * </ul>
 * </p>
 */
public final class Shanks extends SeriesAlgorithm {

	private final int myColN;
	private final double[] myG;
	private final double[] myR;

	private double myTerm;

	public Shanks(final double tolerance, final int maxIters) {
		super(tolerance, maxIters);
		myColN = 20;
		myG = new double[myColN << 1];
		myR = new double[myColN << 1];
	}

	@Override
	public final double next(final double e, final double term) {

		// wait to compute series difference
		if (myIndex == 0) {
			++myIndex;
			myTerm = term;
			return term;
		}

		// initialize table on the first element and difference
		if (myIndex == 1) {
			myR[0] = term - myTerm;
			myG[0] = myTerm;
			++myIndex;
			myTerm = term;
			return term;
		}

		// process the next element
		++myIndex;
		final int idx = Math.min(myIndex - 1, (myColN << 1)) - 1;
		double tmp1 = term - myTerm;
		double tmp2 = 0.0, tmp3 = 0.0;
		double extrap = myTerm;
		for (int i = 0; i <= idx - 1; ++i) {
			if (Math.abs(myR[i]) < TINY) {
				myTerm = term;
				return Double.NaN;
			}
			final double rs = tmp1 / myR[i] - 1.0;
			final double tmp0;
			if (i == 0) {
				tmp0 = rs;
			} else {
				tmp0 = rs * myR[i - 1];
			}
			if ((i & 1) == 0) {
				if (Math.abs(rs) < TINY) {
					myTerm = term;
					return Double.NaN;
				}
				if (i == 0) {
					tmp3 = myG[0] - (extrap - myG[0]) / rs;
					myG[0] = extrap;
					if (myIndex <= 4) {
						extrap = tmp3;
					}
				} else {
					extrap = myG[i - 1] - (myG[i] - myG[i - 1]) / rs;
					myG[i - 1] = myG[i];
					myG[i] = tmp3;
					tmp3 = extrap;
				}
			}
			if (i > 0) {
				myR[i - 1] = tmp2;
			}
			tmp2 = tmp1;
			tmp1 = tmp0;
		}
		myR[idx - 1] = tmp2;
		if (idx < (myColN << 1)) {
			myR[idx] = tmp1;
		}
		myG[idx] = extrap;
		myTerm = term;
		return extrap;
	}

	@Override
	public final String getName() {
		return "Shanks";
	}
}
