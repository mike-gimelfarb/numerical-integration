package algebras;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 * A class of extended precision floating point numbers, implemented according
 * to the IEEE 754R format.
 */
public final class ExtendedPrecision implements OrderedField<BigDecimal> {

    /**
     * A decimal system matching the IEEE 754R Decimal32 format, 7 digits, and a
     * rounding mode of {@link RoundingMode#HALF_EVEN HALF_EVEN}, the IEEE 754R
     * default.
     */
    public static final ExtendedPrecision DECIMAL32 = new ExtendedPrecision(MathContext.DECIMAL32);

    /**
     * A decimal system matching the IEEE 754R Decimal64 format, 16 digits, and a
     * rounding mode of {@link RoundingMode#HALF_EVEN HALF_EVEN}, the IEEE 754R
     * default.
     */
    public static final ExtendedPrecision DECIMAL64 = new ExtendedPrecision(MathContext.DECIMAL64);

    /**
     * A decimal system matching the IEEE 754R Decimal128 format, 34 digits, and a
     * rounding mode of {@link RoundingMode#HALF_EVEN HALF_EVEN}, the IEEE 754R
     * default.
     */
    public static final ExtendedPrecision DECIMAL128 = new ExtendedPrecision(MathContext.DECIMAL128);

    /**
     * A decimal system matching the IEEE 754R Decimal256 format, 71 digits, and a
     * rounding mode of {@link RoundingMode#HALF_EVEN HALF_EVEN}, the IEEE 754R
     * default.
     */
    public static final ExtendedPrecision DECIMAL256 = new ExtendedPrecision(
	    new MathContext(71, RoundingMode.HALF_EVEN));

    protected final MathContext myContext;

    /**
     * Creates a new extended precision instance.
     * 
     * @param context a MathContext object containing the precision and rounding
     *                mode to use for calculations.
     */
    public ExtendedPrecision(final MathContext context) {
	myContext = context;
    }

    @Override
    public final BigDecimal zero() {
	return BigDecimal.ZERO;
    }

    @Override
    public final BigDecimal one() {
	return BigDecimal.ONE;
    }

    @Override
    public final BigDecimal huge() {
	return BigDecimal.valueOf(1e60);
    }

    @Override
    public final BigDecimal negate(final BigDecimal e) {
	return e.negate();
    }

    @Override
    public final BigDecimal invert(final BigDecimal e) {
	return BigDecimal.ONE.divide(e, myContext);
    }

    @Override
    public final BigDecimal add(final BigDecimal a, final BigDecimal b) {
	return a.add(b, myContext);
    }

    @Override
    public final BigDecimal multiply(final BigDecimal a, final BigDecimal b) {
	return a.multiply(b, myContext);
    }

    @Override
    public final double distance(final BigDecimal a, final BigDecimal b) {
	return a.subtract(b, myContext).abs().doubleValue();
    }
}
