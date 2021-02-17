package series.field;

import algebras.OrderedField;

/**
 * Implements the Iterated Theta algorithm for convergence of infinite series as
 * described in [1].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Weniger, Ernst Joachim. "Nonlinear sequence transformations for the
 * acceleration of convergence and the summation of divergent series." arXiv
 * preprint math/0306302 (2003).</li>
 * </ul>
 * </p>
 */
public final class IteratedTheta<T> extends SeriesAlgorithm<T> {

    private final T[] myTab;

    public IteratedTheta(final OrderedField<T> field, final double tolerance, final double underflow,
	    final int maxIters) {
	super(field, tolerance, underflow, maxIters);
	myTab = myField.emptyList(myMaxIters);
    }

    @Override
    public final T next(final T e, final T term) {

	// initial value for J
	myTab[myIndex] = term;
	if (myIndex < 3) {
	    ++myIndex;
	    return term;
	}

	// higher-order J array
	final int lmax = myIndex / 3;
	int m = myIndex;
	for (int l = 1; l <= lmax; ++l) {
	    m -= 3;
	    final T diff0 = myField.subtract(myTab[m + 1], myTab[m]);
	    final T diff1 = myField.subtract(myTab[m + 2], myTab[m + 1]);
	    final T diff2 = myField.subtract(myTab[m + 3], myTab[m + 2]);
	    final T dfsq1 = myField.subtract(diff1, diff0);
	    final T dfsq2 = myField.subtract(diff2, diff1);
	    final T denom = myField.subtract(myField.multiply(diff2, dfsq1), myField.multiply(diff0, dfsq2));

	    // safeguarded against underflow
	    if (myField.magnitude(denom) <= myTiny) {
		myTab[m] = myField.huge();
	    } else {
		final T ratio = myField.divide(dfsq2, denom);
		final T prod = myField.multiply(myField.multiply(diff0, diff1), ratio);
		myTab[m] = myField.subtract(myTab[m + 1], prod);
	    }
	}
	++myIndex;

	// return result
	return myTab[myIndex % 3];
    }

    @Override
    public String getName() {
	return "Iterated Theta";
    }
}
