import java.math.BigDecimal;
import java.util.function.Function;

import algebras.ExtendedPrecision;
import series.field.*;

public class series_test_extended {

	public static double factorial(long n) {
		if (n <= 1) {
			return 1.;
		} else {
			return n * factorial(n - 1);
		}
	}

	public static void main(final String[] args) {

		ExtendedPrecision prec = ExtendedPrecision.DECIMAL256;
		Function<Long, BigDecimal> func = k -> prec.invert(prec.multiply(BigDecimal.valueOf(k), BigDecimal.valueOf(k)));
//
//		SeriesAlgorithm alg = new Adaptive(1e-6, 10000, 15, false);
//		Function<Long, Double> alt = alg.alternatingSeries(func, 1L);

		SeriesAlgorithm<BigDecimal> alg2 = new BrezinskiTheta<>(prec, 1e-40, 1e-200, 190000);
		BigDecimal lim = alg2.limit(func, 1L, true);
		System.out.println(lim);
	}
}
