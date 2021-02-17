package integral;

import java.util.Stack;
import java.util.function.Function;

/**
 * A globally adaptive numerical integrator based on Simpson's rule [1]. The
 * Lyness-Richardson criterion is used to decide when to stop subdividing an
 * interval [2].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Guy F. Kuncir. 1962. Algorithm 103: Simpson's rule integrator.
 * Commun. ACM 5, 6 (June 1962), 347.
 * DOI:https://doi.org/10.1145/367766.368179</li>
 * <li>[2] Lyness, James N. "Notes on the adaptive Simpson quadrature routine."
 * Journal of the ACM (JACM) 16.3 (1969): 483-495.</li>
 * </ul>
 * </p>
 */
public final class Simpson extends Quadrature {

    private final int myMaxEvals;

    /**
     * Creates a new instance of the adaptive Simpson integrator.
     * 
     * @param tolerance the smallest acceptable change in integral estimates in
     *                  consecutive iterations that indicates the algorithm has
     *                  converged
     * @param maxEvals  the maximum number of function evaluations
     */
    public Simpson(final double tolerance, final int maxEvals) {
	super(tolerance);
	myMaxEvals = maxEvals;
    }

    public Simpson(final double tolerance) {
	this(tolerance, 9999);
    }

    @Override
    final double properIntegral(final Function<? super Double, Double> f, final double a, final double b) {
	return adaptiveSimpson(f, a, b);
    }

    private double adaptiveSimpson(final Function<? super Double, Double> func, final double a, final double b) {

	// INITIALIZE ITERATION VARIABLES
	final Stack<double[]> stack = new Stack<>();

	// ADD THE INITIAL INTERVAL
	final double m0 = 0.5 * (a + b);
	final double fa = func.apply(a);
	final double fm = func.apply(m0);
	final double fb = func.apply(b);
	double[] interval = { a, m0, b, fa, fm, fb };
	stack.push(interval);
	int fev = 3;
	int memsize = 1;
	double est = 0.0;

	// MAIN LOOP OF GLOBALLY ADAPTIVE SIMPSON
	while (!stack.isEmpty()) {
	    memsize = Math.max(memsize, stack.size());

	    // RETRIEVE AND REMOVE A SUBINTERVAL FROM THE QUEUE
	    interval = stack.pop();
	    final double a1 = interval[0];
	    final double m1 = interval[1];
	    final double b1 = interval[2];
	    final double fa1 = interval[3];
	    final double fm1 = interval[4];
	    final double fb1 = interval[5];

	    // ESTIMATE ERROR USING LYNESS' RICHARDSON METHOD
	    final double wid = b1 - a1;
	    final double lm1 = 0.5 * (a1 + m1);
	    final double rm1 = 0.5 * (m1 + b1);
	    final double flm1 = func.apply(lm1);
	    final double frm1 = func.apply(rm1);
	    final double est1 = (fa1 + 4.0 * fm1 + fb1) * (wid / 6.0);
	    final double est2 = (fa1 + 4.0 * flm1 + fm1) * (wid / 12.0) + (fm1 + 4.0 * frm1 + fb1) * (wid / 12.0);
	    fev += 2;

	    // IF THE ERROR IS WITHIN TOLERANCE, UPDATE ESTIMATE, ELSE DEQUEUE
	    final double err = Math.abs(est1 - est2) / 15.0;
	    if (err <= wid * myTol) {
		est += est2 + (est2 - est1) / 15.0;
		continue;
	    }

	    // CHECK FOR FAILURE OF CONVERGENCE
	    if (fev >= myMaxEvals || !Double.isFinite(est)) {
		myFEvals += fev;
		return Double.NaN;
	    }

	    // ADD THE LEFT HALF-INTERVAL
	    final double[] left = { a1, lm1, m1, fa1, flm1, fm1 };
	    stack.push(left);

	    // ADD THE RIGHT HALF-INTERVAL
	    final double[] right = { m1, rm1, b1, fm1, frm1, fb1 };
	    stack.push(right);
	}

	// RETURN THE FINAL ESTIMATE
	myFEvals += fev;
	return est;
    }
}
