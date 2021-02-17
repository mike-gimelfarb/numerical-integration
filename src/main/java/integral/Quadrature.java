package integral;

import java.util.Iterator;
import java.util.function.Function;

import series.dp.SeriesAlgorithm;

/**
 * An abstract class for the base of all numerical integrators for 1d functions.
 */
public abstract class Quadrature {

    protected static final int MAX_ITERS = (int) 1e7;

    protected final double myTol;
    protected int myFEvals;

    /**
     * Creates a new instance of the current numerical integrator.
     * 
     * @param tolerance the smallest acceptable change in integral estimates in
     *                  consecutive iterations that indicates the algorithm has
     *                  converged
     */
    public Quadrature(final double tolerance) {
	myTol = tolerance;
	myFEvals = 0;
    }

    abstract double properIntegral(Function<? super Double, Double> f, double a, double b);

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
	    return -integrate(f, b, a);
	}

	// finite integral (a, b)
	if (Double.isFinite(a) && Double.isFinite(b)) {
	    return properIntegral(f, a, b);
	}

	// improper integral is not handled by default
	throw new IllegalArgumentException("Quadrature subroutine" + " " + getClass().getSimpleName() + " "
		+ "does not handle improper integrals. " + "Use Gauss Kronrod subroutine instead "
		+ "or integrate piecewise.");
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
     * Computes an indefinite integral of an oscillating function (a function whose
     * sign changes periodically) from a given starting point to infinity. This
     * method first splits the integration region into disjoint regions. It then
     * defines a sequence by estimating the definite integrals over each disjoint
     * region, and accelerates the sequence of partial sums to estimate the
     * indefinite integral.
     * 
     * @param f               the function to integrate
     * @param a               the left endpoint of the integration region
     * @param splitPoints     the points at which to split the integration region
     *                        when defining definite integrals
     * @param maxIntegrations the maximum number of integrations that can be
     *                        performed to estimate the indefinite integral
     * @return an estimate of the indefinite integral
     */
    public double integratePiecewise(final Function<? super Double, Double> f, final double a,
	    final int maxIntegrations, final Iterable<Double> splitPoints, final SeriesAlgorithm accelerator) {

	// cycle through the roots until we find the first in the integration region
	final Iterator<Double> it = splitPoints.iterator();
	double x0;
	while ((x0 = it.next()) <= a)
	    ;

	// get the first root and integral estimate in [a, x(0)]
	final double a1 = x0;
	final double ix0 = properIntegral(f, a, a1);

	// construct a sequence of integrals over the intervals between successive roots
	final Iterable<Double> sequence = () -> new Iterator<>() {

	    double left = 0.0;
	    double right = a1;
	    int i = 1;

	    @Override
	    public boolean hasNext() {
		return i < maxIntegrations && it.hasNext();
	    }

	    @Override
	    public Double next() {
		left = right;
		right = it.next();
		++i;
		return properIntegral(f, left, right);
	    }
	};

	// accelerate the sequence
	return ix0 + accelerator.limit(sequence, true);
    }
}
