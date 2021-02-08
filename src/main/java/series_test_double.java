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

		Function<Long, Double> func = k -> Math.pow(-1, k) / k;
//
//		SeriesAlgorithm alg = new Adaptive(1e-6, 10000, 15, false);
//		Function<Long, Double> alt = alg.alternatingSeries(func, 1L);

		SeriesAlgorithm alg2 = new WynnRho(1e-7, 10000);
		double lim = alg2.limit(func, 1L, true);
		System.out.println(lim);
		System.out.println(alg2.countEvaluations());
	}
}
