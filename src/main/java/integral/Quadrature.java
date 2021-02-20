package integral;

import java.util.Iterator;
import java.util.function.Function;

/**
 * An abstract class for the base of all numerical integrators for 1d functions.
 */
public abstract class Quadrature {

    protected final double myTol;
    protected final int myMaxEvals;

    protected int myFEvals;

    /**
     * Creates a new instance of the current numerical integrator.
     * 
     * @param tolerance      the smallest acceptable change in integral estimates in
     *                       consecutive iterations that indicates the algorithm has
     *                       converged
     * @param maxEvaluations the maximum number of evaluations of each function
     *                       permitted
     */
    public Quadrature(final double tolerance, final int maxEvaluations) {
	myTol = tolerance;
	myMaxEvals = maxEvaluations;
	myFEvals = 0;
    }

    abstract double properIntegral(Function<? super Double, Double> f, double a, double b);

    /**
     * Returns a string representation of the current algorithm.
     * 
     * @return a string representation of the current algorithm
     */
    public abstract String getName();

    /**
     * Sets the variable that keeps track of the number of function evaluations to
     * zero.
     */
    public final void resetCounter() {
	myFEvals = 0;
    }

    /**
     * Returns the number of evaluations of the function made thus far.
     * 
     * @return an integer representing the number of evaluations of the function
     *         made thus far.
     */
    public final int getEvaluations() {
	return myFEvals;
    }

    /**
     * Returns the tolerance of this instance.
     * 
     * @return the smallest acceptable change in integral estimates in consecutive
     *         iterations that indicates the algorithm has converged
     */
    public final double getTolerance() {
	return myTol;
    }

    /**
     * Estimates a definite integral, or indefinite integral in some limited cases.
     * All algorithms support definite integrals. The only algorithm that currently
     * supports indefinite integrals is Gauss Kronrod.
     * 
     * @param f the function to integrate
     * @param a the left endpoint of the integration interval
     * @param b the right endpoint of the integration interval
     * @return an estimate of the definite or indefinite integral
     */
    public double integrate(final Function<? super Double, Double> f, final double a, final double b) {

	// empty integral (a, a)
	if (a == b) {
	    return 0.0;
	}

	// opposite bounds
	if (a > b) {
	    return integrate(f, b, a);
	}

	// finite integral (a, b)
	if (Double.isFinite(a) && Double.isFinite(b)) {
	    return properIntegral(f, a, b);
	}

	// improper integral is not handled by default
	throw new IllegalArgumentException("The quadrature subroutine" + " " + getClass().getSimpleName() + " "
		+ "does not handle improper integrals. Try using a transformation or integrate piecewise.");
    }

    /**
     * Estimates a definite integral of a function by first applying a change of
     * variable to the function, suitable for estimating indefinite integrals. Given
     * a function to integrate, f, and a differentiable function t, this method
     * integrates f(t(x)) * t'(x).
     * 
     * @param f  the function to integrate
     * @param t  the change of variable
     * @param dt the derivative of the change of variable function t
     * @param a  the left endpoint of the integration interval
     * @param b  the right endpoint of the integration interval
     * @return an estimate of the definite integral
     */
    public double integrate(final Function<? super Double, Double> f, final Function<? super Double, Double> t,
	    final Function<? super Double, Double> dt, final double a, final double b) {

	// construct the integrand f(t(x)) t(x) dt(x)
	final Function<Double, Double> func = (x) -> f.apply(t.apply(x)) * dt.apply(x);

	// dispatch to integration routine
	return integrate(func, a, b);
    }

    /**
     * This method splits a potentially indefinite integration interval into
     * disjoint definite intervals. It returns a sequence of integral estimates on
     * each interval. This is very useful for evaluating highly oscillating
     * integrals on infinite intervals, called Longman's method [1].
     * 
     * <p>
     * References:
     * <ul>
     * <li>[1] Longman, I. (1956). Note on a method for computing infinite integrals
     * of oscillatory functions. Mathematical Proceedings of the Cambridge
     * Philosophical Society, 52(4), 764-768. doi:10.1017/S030500410003187X</li>
     * </ul>
     * </p>
     * 
     * @param f           the function to integrate
     * @param splitPoints the points at which to split the integration region when
     *                    defining definite integrals
     * @return an Iterable representing the sequence of definite integrals
     */
    public Iterable<Double> integratePiecewise(final Function<? super Double, Double> f,
	    final Iterable<Double> splitPoints) {
	return () -> new Iterator<>() {

	    private final Iterator<Double> it = splitPoints.iterator();
	    private double left = 0.0;
	    private double right = it.hasNext() ? it.next() : Double.NaN;

	    @Override
	    public final boolean hasNext() {
		return it.hasNext();
	    }

	    @Override
	    public final Double next() {
		left = right;
		right = it.next();
		return integrate(f, left, right);
	    }
	};
    }
}
