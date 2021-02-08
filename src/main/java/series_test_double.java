import java.math.BigDecimal;
import java.util.function.Function;

import series.dp.*;

public class series_test_double {

	public static double factorial(long n) {
		if (n <= 1) {
			return 1.;
		} else {
			return n * factorial(n - 1);
		}
	}

	public static void main(final String[] args) {

		Function<Long, Double> func = k -> 1. / Math.pow(k, 2);
//
//		SeriesAlgorithm alg = new Adaptive(1e-6, 10000, 15, false);
//		Function<Long, Double> alt = alg.alternatingSeries(func, 1L);

		SeriesAlgorithm alg2 = new IteratedTheta(1e-6, 10000);
		for (int i = 1; i <= 1000; ++i) {
			double lim = alg2.limit(func, 1L, true);
			System.out.println(lim);
		}
	}
}
