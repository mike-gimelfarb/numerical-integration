package algebras;

import java.lang.reflect.Array;

/**
 * An abstract ordered field with a metric.
 * 
 * @param <T>the type of the objects in the field
 */
public interface OrderedField<T> {

    /**
     * Returns the zero element of the field.
     * 
     * @return the zero element
     */
    public T zero();

    /**
     * Returns the identity element of the field.
     * 
     * @return the identity element
     */
    public T one();

    /**
     * Returns a number that is considered very large in the field. Used to
     * safeguard for underflow in calculations.
     * 
     * @return a large element of the field
     */
    public T huge();

    /**
     * Returns the additive inverse of an element of the field.
     * 
     * @param e the element to invert
     * @return the additive inverse
     */
    public T negate(T e);

    /**
     * Returns the multiplicative inverse of an element in the field.
     * 
     * @param e the element to invert
     * @return the multiplicative inverse
     */
    public T invert(T e);

    /**
     * Performs an addition of two elements in the field according to the additive
     * operator.
     * 
     * @param a the first element
     * @param b the second element
     * @return the sum of the two elements
     */
    public T add(T a, T b);

    /**
     * Performs a multiplication of two elements in the field according to the
     * multiplicative operator.
     * 
     * @param a the first element
     * @param b the second element
     * @return the product of the two elements
     */
    public T multiply(T a, T b);

    /**
     * Computes the distance between two elements in the field according to an
     * induced metric.
     * 
     * @param a the first element
     * @param b the second element
     * @return the distance between the two elements
     */
    public double distance(T a, T b);

    /**
     * Computes the magnitude of an element in the field according to the induced
     * metric. This is defined as d(e, 0), where 0 is the zero element of the field
     * and d is the distance metric.
     * 
     * @param e the element
     * @return the magnitude of the element
     */
    default public double magnitude(final T e) {
	return distance(e, zero());
    }

    /**
     * Computes the additive difference between two elements in the field. This is
     * defined as the sum of the first element and the additive inverse of the
     * second element.
     * 
     * @param a the first element
     * @param b the second element
     * @return the difference between the two elements
     */
    default public T subtract(final T a, final T b) {
	return add(a, negate(b));
    }

    /**
     * Computes the quotient between two elements in the field. This is defined as
     * the product of the first element and the multiplicative inverse of the second
     * element.
     * 
     * @param a the first element
     * @param b the second element
     * @return the quotient of the two elements
     */
    default public T divide(final T a, final T b) {
	return multiply(a, invert(b));
    }

    /**
     * Returns an empty array of the specified length that can store items in this
     * field. The array initially contains null elements.
     * 
     * @param length the desired length of the array
     * @return an empty array of elements of type T
     */
    @SuppressWarnings("unchecked")
    default public T[] emptyList(final int length) {
	return (T[]) Array.newInstance(zero().getClass(), length);
    }
}
