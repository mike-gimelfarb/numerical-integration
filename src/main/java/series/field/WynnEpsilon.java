package series.field;

import algebras.OrderedField;

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
public final class WynnEpsilon<T> extends SeriesAlgorithm<T> {

    private final T[] myTab;

    public WynnEpsilon(final OrderedField<T> field, final double tolerance, final double underflow,
	    final int maxIters) {
	super(field, tolerance, underflow, maxIters);
	myTab = myField.emptyList(maxIters + 1);
    }

    @Override
    public final T next(final T e, final T term) {

	// initialization
	myTab[myIndex] = term;
	if (myIndex == 0) {
	    ++myIndex;
	    return term;
	}

	// main loop of Wynn Epsilon algorithm
	T aux = myField.zero();
	for (int j = myIndex; j >= 1; --j) {

	    // compute the next epsilon
	    final T temp = aux;
	    aux = myTab[j - 1];
	    T diff = myField.subtract(myTab[j], aux);

	    // correct denominators for underflow
	    if (myField.magnitude(diff) <= myTiny) {
		diff = myField.huge();
	    } else {
		diff = myField.add(temp, myField.invert(diff));
	    }
	    myTab[j - 1] = diff;
	}

	// prepare result
	final T result = myTab[myIndex & 1];
	++myIndex;
	return result;
    }

    @Override
    public String getName() {
	return "Wynn Epsilon";
    }
}
