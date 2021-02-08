package series.dp;

import java.util.ArrayList;
import java.util.List;

import series.dp.Levin.RemainderSequence;

/**
 * An adaptive algorithm for evaluating the limits of sequences and series of
 * real values. This algorithms works by initializing a number of convergence
 * acceleration algorithms with different properties, and running them in
 * parallel. The estimates of the sequence's or corresponding series limit are
 * extracted from the instance that has converged first to a stable limit to the
 * desired tolerance.
 */
public final class Adaptive extends SeriesAlgorithm {

	private final boolean myPrint;
	private final int myConvergeIters;
	private final List<SeriesAlgorithm> myMethods;
	private final List<Boolean> myForAlternating;

	/**
	 * Creates a new instance of an adaptive convergence acceleration algorithm.
	 * 
	 * @param tolerance                the smallest acceptable change in series
	 *                                 evaluation in consecutive iterations that
	 *                                 indicates the algorithm has converged
	 * @param maxIters                 the maximum number of sequence terms to
	 *                                 evaluate before giving up
	 * @param iterationsForConvergence the number of iterations used to assess
	 *                                 whether or not a particular instance of
	 *                                 convergence acceleration algorithm used by
	 *                                 this algorithm has converged
	 * @param printProgress            whether or not to print the progress of the
	 *                                 current algorithm, e.g. the estimate of a
	 *                                 sequence's or series's limit at each
	 *                                 iteration according to each algorithm
	 *                                 instance
	 */
	public Adaptive(final double tolerance, final int maxIters, final int iterationsForConvergence,
			final boolean printProgress) {
		super(tolerance, maxIters);
		myPrint = printProgress;
		myConvergeIters = iterationsForConvergence;
		myMethods = new ArrayList<>();
		myForAlternating = new ArrayList<>();

		// for alternating series
		myMethods.add(new Cohen(myTol, myMaxIters));
		myForAlternating.add(true);

		// for series of positive or negative terms
		// for accelerating linear and logarithmic convergence
		myMethods.add(new Levin(myTol, myMaxIters, RemainderSequence.T));
		myMethods.add(new Levin(myTol, myMaxIters, RemainderSequence.U));
		myMethods.add(new Levin(myTol, myMaxIters, RemainderSequence.V));
		myMethods.add(new BrezinskiTheta(myTol, myMaxIters));
		myMethods.add(new IteratedTheta(myTol, myMaxIters));
		myMethods.add(new WynnRho(myTol, myMaxIters));
		myForAlternating.add(false);
		myForAlternating.add(false);
		myForAlternating.add(false);
		myForAlternating.add(false);
		myForAlternating.add(false);
		myForAlternating.add(false);

		// for series with irregular signs
		myMethods.add(new WynnEpsilon(myTol, myMaxIters));
		myForAlternating.add(false);
	}

	/**
	 * Creates a new instance of an adaptive convergence acceleration algorithm.
	 * 
	 * @param tolerance                the smallest acceptable change in series
	 *                                 evaluation in consecutive iterations that
	 *                                 indicates the algorithm has converged
	 * @param maxIters                 the maximum number of sequence terms to
	 *                                 evaluate before giving up
	 * @param iterationsForConvergence the number of iterations used to assess
	 *                                 whether or not a particular instance of
	 *                                 convergence acceleration algorithm used by
	 *                                 this algorithm has converged
	 */
	public Adaptive(final double tolerance, final int maxIters, final int iterationsForConvergence) {
		this(tolerance, maxIters, iterationsForConvergence, false);
	}

	@Override
	public final double next(final double e, final double term) {
		// this method should not be applied sequentially
		return Double.NaN;
	}

	@Override
	public final double limit(final Iterable<? extends Double> seq, final boolean series) {

		// reset the storage counter and initialize variables for result
		myIndex = 0;
		for (final SeriesAlgorithm method : myMethods) {
			method.myIndex = 0;
		}

		// initialize temporary variables
		double term = 0.0;
		final int mcount = myMethods.size();
		final double[] oldests = new double[mcount];
		final double[] ests = new double[mcount];
		final boolean[] except = new boolean[mcount];
		final int[] converges = new int[mcount];

		// initialize variables for tracking sign
		double oldsignum = 0.0;
		double globalsignum = 0.0;
		boolean alternates = true;

		// print header
		final int printPrec = 18;
		if (myPrint) {
			String line = "";
			line += pad("Index", printPrec + 5) + "\t";
			for (final SeriesAlgorithm method : myMethods) {
				line += pad(method.getName(), printPrec + 5) + "\t";
			}
			System.out.println(line);
		}

		// main loop
		for (final Double e : seq) {
			++myFEvals;
			if (e == null || !Double.isFinite(e)) {
				if (myPrint) {
					System.out.println("Aborting series acceleration" + " " + "due to NaN or null" + " "
							+ "term at iteration" + " " + myIndex);
				}
				break;
			}

			// get the next term of the sequence
			if (series) {
				term += e;
			} else {
				term = e;
			}

			// track the sign of the sequence
			final boolean oldalternates = alternates;
			final double signum = Math.signum(e);
			if (myIndex == 0) {
				globalsignum = signum;
			} else {
				if (signum != globalsignum) {
					globalsignum = 0.0;
				}
				if (alternates && signum == oldsignum) {
					alternates = false;
				}
			}
			oldsignum = signum;

			// exclude alternating series methods
			if (alternates != oldalternates) {
				for (int m = 0; m < mcount; ++m) {
					final SeriesAlgorithm method = myMethods.get(m);
					if (myForAlternating.get(m)) {
						except[m] = true;
						if (myPrint) {
							System.out.println("Disabling method" + " " + method.getName() + " "
									+ "since series is nonalternating" + " " + "at iteration" + " " + myIndex);
						}
					}
				}
			}

			// this is the estimation and checking stage
			// note that we simply iterate over all methods and sequentially
			// estimate the next term and error from each, and check convergence
			for (int m = 0; m < mcount; ++m) {

				// get the method at index m
				final SeriesAlgorithm method = myMethods.get(m);
				if (except[m]) {
					continue;
				}

				// estimate the next terms and errors
				oldests[m] = ests[m];
				ests[m] = method.next(e, term);
				final double error = Math.abs(oldests[m] - ests[m]);

				// convergence test
				if (Math.abs(ests[m]) < HUGE && error <= myTol) {
					++converges[m];

					// has the current method converged?
					if (converges[m] >= myConvergeIters) {
						if (myPrint) {
							System.out.println("Converged after" + " " + myIndex + " " + "iterations with method" + " "
									+ method.getName());
						}
						return ests[m];
					}
				} else {
					converges[m] = 0;
				}

				// failed to converge in maximum number of iterations
				if (myIndex >= myMaxIters) {
					return Double.NaN;
				}

				// if a method produces an invalid estimate it's excluded
				if (myIndex > myConvergeIters && !Double.isFinite(ests[m])) {
					if (myPrint) {
						System.out.println("Disabling method" + " " + method.getName() + " " + "due to instability"
								+ " " + "at iteration" + " " + myIndex);
					}
					except[m] = true;
				}
			}
			++myIndex;

			// print progress
			if (myPrint) {
				String line = pad(myIndex + "", printPrec) + "\t";
				for (int m = 0; m < mcount; ++m) {
					if (except[m]) {
						line += pad("-", printPrec) + "\t";
					} else {
						line += pad(Double.toString(ests[m]), printPrec) + "\t";
					}
				}
				System.out.println(line);
			}
		}

		// did not achieve the desired error
		return Double.NaN;
	}

	@Override
	public String getName() {
		return "Adaptive";
	}

	private static String pad(final String str, final int len) {
		if (str.length() >= len) {
			return str;
		}
		String result = str;
		for (int i = str.length() + 1; i <= len; ++i) {
			result += " ";
		}
		return result;
	}
}
