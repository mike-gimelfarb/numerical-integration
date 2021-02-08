package series.field;

import algebras.OrderedField;

/**
 * Implements Brezinski's Theta algorithm for convergence of infinite series
 * based on [1] and as described in [2].
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Brezinski, Claude. Accélération de la convergence en analyse
 * numérique. Vol. 584. Springer, 2006.</li>
 * <li>[2] Weniger, Ernst Joachim. "Nonlinear sequence transformations for the
 * acceleration of convergence and the summation of divergent series." arXiv
 * preprint math/0306302 (2003).</li>
 * </ul>
 * </p>
 */
public final class BrezinskiTheta<T> extends SeriesAlgorithm<T> {

	private T[] myTabA, myTabB;

	public BrezinskiTheta(final OrderedField<T> field, final double tolerance, final double underflow,
			final int maxIters) {
		super(field, tolerance, underflow, maxIters);
		myTabA = myField.emptyList(maxIters);
		myTabB = myField.emptyList(maxIters);
	}

	@Override
	public final T next(final T e, final T term) {

		// first entry of table
		if (myIndex == 0) {
			myTabA[0] = term;
			++myIndex;
			return term;
		}

		// swapping the A and B arrays
		final T[] a, b;
		if ((myIndex & 1) == 0) {
			a = myTabA;
			b = myTabB;
		} else {
			a = myTabB;
			b = myTabA;
		}

		// the case n >= 1
		final int jmax = ((myIndex << 1) + 1) / 3;
		T aux1 = a[0];
		T aux2 = myField.zero();
		a[0] = term;
		for (int j = 1; j <= jmax; ++j) {
			final T aux3 = aux2;
			aux2 = aux1;
			if (j < jmax) {
				aux1 = a[j];
			}
			if ((j & 1) == 0) {
				final T twob = myField.add(b[j - 1], b[j - 1]);
				final T denom = myField.add(myField.subtract(a[j - 1], twob), aux2);

				// correct for underflow
				if (myField.magnitude(denom) <= myTiny) {
					a[j] = myField.huge();
				} else {
					final T diff1 = myField.subtract(b[j - 2], aux3);
					final T diff2 = myField.subtract(a[j - 1], b[j - 1]);
					final T delta = myField.divide(myField.multiply(diff1, diff2), denom);
					a[j] = myField.add(aux3, delta);
				}
			} else {
				final T diff = myField.subtract(a[j - 1], b[j - 1]);

				// correct for underflow
				if (myField.magnitude(diff) <= myTiny) {
					a[j] = myField.huge();
				} else {
					a[j] = myField.add(aux3, myField.invert(diff));
				}
			}
		}
		++myIndex;

		// return result
		if ((jmax & 1) == 0) {
			return a[jmax];
		} else {
			return a[jmax - 1];
		}
	}

	@Override
	public String getName() {
		return "Brezinski Theta";
	}
}
