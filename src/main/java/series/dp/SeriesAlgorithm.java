package series.dp;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.function.Function;

/**
 * The base class of all convergence acceleration algorithms. Provides
 * functionality for evaluating infinite series and limits of sequences of real
 * values.
 */
public abstract class SeriesAlgorithm {

	protected static final double HUGE = 1e60;
	protected static final double TINY = 1e-60;

	protected final double myTol;
	protected final int myMaxIters;
	protected int myFEvals;
	protected int myIndex;

	/**
	 * Creates a new instance of the current convergence acceleration algorithm.
	 * 
	 * @param tolerance the smallest acceptable change in series evaluation in
	 *                  consecutive iterations that indicates the algorithm has
	 *                  converged
	 * @param maxIters  the maximum number of sequence terms to evaluate before
	 *                  giving up
	 */
	public SeriesAlgorithm(final double tolerance, final int maxIters) {
		myTol = tolerance;
		myMaxIters = maxIters;
		myFEvals = 0;
	}

	/**
	 * Update the current algorithm after observing the specified element and
	 * partial sum.
	 * 
	 * @param e    the next element of the sequence
	 * @param term the partial sum of all observed elements of the sequence thus far
	 * @return the next estimate of the limit of the sequence
	 */
	public abstract double next(double e, double term);

	/**
	 * Returns a string representation of the current algorithm.
	 * 
	 * @return a string representation of the current algorithm
	 */
	public abstract String getName();

	/**
	 * Sets the variable that keeps track of the number of sequence evaluations to
	 * zero.
	 */
	public final void resetCounter() {
		myFEvals = 0;
	}

	/**
	 * Returns the number of evaluations of the sequence made thus far.
	 * 
	 * @return an integer representing the number of evaluations of the sequence
	 *         made thus far.
	 */
	public final int countEvaluations() {
		return myFEvals;
	}

	/**
	 * Given a sequence represented as an Iterable object, numerically evaluates the
	 * limit of the sequence or the corresponding series.
	 * 
	 * @param seq    an Iterable object of type Double whose limit to evaluate
	 * @param series a boolean variable indicating whether to evaluate the limit of
	 *               the series or the sequence
	 * @return a Double that approximates the limit of the sequence or corresponding
	 *         series. If the limit cannot be determined, returns NaN
	 */
	public double limit(final Iterable<? extends Double> seq, final boolean series) {
		myIndex = 0;
		double term = 0.0;
		double est = Double.NaN;

		for (final double e : seq) {

			// get the next term of the sequence
			if (series) {
				term += e;
			} else {
				term = e;
			}
			++myFEvals;

			// estimate the next term and error
			final double oldest = est;
			est = next(e, term);
			if (myIndex >= 4) {
				if (est != est) {
					break;
				}
				final double error = Math.abs(oldest - est);
				if (error <= myTol) {
					return est;
				}
				if (myIndex >= myMaxIters || !Double.isFinite(error)) {
					break;
				}
			}
		}

		// did not achieve the desired error
		return Double.NaN;
	}

	/**
	 * Given a sequence represented as a function and a starting index, numerically
	 * evaluates the limit of the sequence or the corresponding series.
	 * 
	 * @param func   a Function representing the sequence
	 * @param start  an long representing the starting index at which the sequence
	 *               begins
	 * @param series a boolean variable indicating whether to evaluate the limit of
	 *               the series or the sequence
	 * @return a Double that approximates the limit of the sequence or corresponding
	 *         series. If the limit cannot be determined, returns NaN
	 */
	public double limit(final Function<? super Long, ? extends Double> func, final long start, final boolean series) {

		// convert function to iterable
		final Iterable<? extends Double> seq = () -> new Iterator<Double>() {

			private long k = start;
			private final Function<? super Long, ? extends Double> myFunc = func;

			@Override
			public boolean hasNext() {
				return k < Long.MAX_VALUE - 1;
			}

			@Override
			public Double next() {
				return myFunc.apply(k++);
			}
		};

		// dispatch
		return limit(seq, series);
	}

	/**
	 * Given a sequence represented as a function and a starting index, numerically
	 * evaluates the limit of the sequence or the corresponding series.
	 * 
	 * @param func   a Function representing the sequence
	 * @param start  an integer representing the starting index at which the
	 *               sequence begins
	 * @param series a boolean variable indicating whether to evaluate the limit of
	 *               the series or the sequence
	 * @return a Double that approximates the limit of the sequence or corresponding
	 *         series. If the limit cannot be determined, returns NaN
	 */
	public double limit(final Function<? super Integer, ? extends Double> func, final int start, final boolean series) {

		// convert function to iterable
		final Iterable<? extends Double> seq = () -> new Iterator<Double>() {

			private int k = start;
			private final Function<? super Integer, ? extends Double> myFunc = func;

			@Override
			public boolean hasNext() {
				return k < Integer.MAX_VALUE - 1;
			}

			@Override
			public Double next() {
				return myFunc.apply(k++);
			}
		};

		// dispatch
		return limit(seq, series);
	}

	/**
	 * Given a sequence represented as a function and a starting index, numerically
	 * evaluates the limit of the sequence or the corresponding series.
	 * 
	 * @param func   a Function representing the sequence
	 * @param start  a BigInteger representing the starting index at which the
	 *               sequence begins
	 * @param series a boolean variable indicating whether to evaluate the limit of
	 *               the series or the sequence
	 * @return a Double that approximates the limit of the sequence or corresponding
	 *         series. If the limit cannot be determined, returns NaN
	 */
	public double limit(final Function<? super BigInteger, ? extends Double> func, final BigInteger start,
			final boolean series) {

		// convert function to iterable
		final Iterable<? extends Double> seq = () -> new Iterator<Double>() {

			private BigInteger k = start;
			private final Function<? super BigInteger, ? extends Double> myFunc = func;

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public Double next() {
				final Double result = myFunc.apply(k);
				k = k.add(BigInteger.ONE);
				return result;
			}
		};

		// dispatch
		return limit(seq, series);
	}

	/**
	 * The Van Wijngaarden transformation converts a convergent sequence of positive
	 * terms into a sequence of alternating terms that has the same limit.
	 * <p>
	 * Formally, given a sequence of terms a(1), a(2)... converts this into an
	 * alternating sequence where the kth term is
	 * <p>
	 * a(k) + 2a(2k) + 4a(4k) + ...
	 * </p>
	 * This method uses the current algorithm's implementation to evaluate the limit
	 * of the above expression. The resulting series are evalated lazily and
	 * represented internally as a function.
	 * </p>
	 * <p>
	 * This trick is incredibly useful for accelerating convergence of
	 * slowly-increasing sequences of positive terms, because alternating series are
	 * generally easier to accelerate.
	 * </p>
	 * <p>
	 * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling,
	 * and Brian P. Flannery. 2007. Numerical Recipes 3rd Edition: The Art of
	 * Scientific Computing (3rd. ed.). Cambridge University Press, USA.
	 * </p>
	 * 
	 * @param seq   a sequence of non-negative terms represented as a Function
	 * @param start a long representing the starting index at which the sequence
	 *              begins
	 * @return a Function representing the resulting alternating sequence obtained
	 *         with Van Wijngaarden's transformation
	 */
	public Function<Long, Double> alternatingSeries(final Function<? super Long, ? extends Double> seq,
			final long start) {
		return (k) -> {

			// create the condensation sequence with kth term 2^k a(2^k)
			final Iterable<Double> seqterm = () -> new Iterator<Double>() {

				private long i = k;
				private double coeff = 1.0;

				@Override
				public final boolean hasNext() {
					return coeff <= (Long.MAX_VALUE >> 1L);
				}

				@Override
				public final Double next() {
					final double term = seq.apply(i + start - 1L);
					final double result = coeff * term;
					i <<= 1L;
					coeff *= 2.0;
					++myFEvals;
					return result;
				}
			};

			// determine the next term as the infinite series of this sequence
			final double term = limit(seqterm, true);
			if (Double.isNaN(term) || ((k - 1L) & 1L) == 0L) {
				return term;
			} else {
				return -term;
			}
		};
	}
}
