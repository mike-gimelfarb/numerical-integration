package test;

import java.util.function.Function;

import series.dp.*;
import series.dp.Levin.RemainderSequence;
import utils.Constants;

public final class SeriesTest {

	static final class InfiniteSeries {

		final Function<Long, Double> mySeries;
		final double myLimit;

		public InfiniteSeries(final Function<Long, Double> series, final double limit) {
			mySeries = series;
			myLimit = limit;
		}
	}

	public static final InfiniteSeries[] ALTERNATING_SERIES = {
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) / (2 * k - 1.0), Math.PI / 4.0),
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) / Math.pow(k, 0.25), 0.5544873859140731),
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) / Math.pow(k, 0.5), 0.60489864342163037),
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) / k, Math.log(2.0)),
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) / Math.pow(k, 1.5), 0.76514702462540794),
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) * Math.log(k + 1) / (k + 1), 0.1598689037),
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) * Math.log(k + 1) / Math.sqrt(k + 1), 0.1932888316),
			new InfiniteSeries(k -> Math.pow(-1.0, k - 1) / Math.log(k + 1), 0.9242998972229) };

	public static final InfiniteSeries[] IRREGULAR_SERIES = {
			new InfiniteSeries(k -> Math.sin(k) / k, 1.0707963267948966),
			new InfiniteSeries(k -> Math.sin(5.0 * k) / k, -0.92920367320510338),
			new InfiniteSeries(k -> Math.sin(k) / k / Math.pow(2.0, k), 0.522937836179332),
			new InfiniteSeries(k -> Math.sin(5.0 * k) / k / Math.pow(2.0, k), -0.509500940834413) };

	public static final InfiniteSeries[] LINEAR_SERIES = {
			new InfiniteSeries(k -> Math.pow(0.8, k) / k, 1.609437912434100),
			new InfiniteSeries(k -> Math.pow(0.995, k - 1), 200.),
			new InfiniteSeries(k -> Math.pow(0.4, k - 1) + Math.pow(0.8, k - 1), 6.666666666666666) };

	public static final InfiniteSeries[] LOGARITHMIC_SERIES = {
			new InfiniteSeries(k -> 1.0 / Math.pow(k, 1.25), 4.59511182584294338),
			new InfiniteSeries(k -> 1.0 / Math.pow(k, 1.5), 2.612375348685488343),
			new InfiniteSeries(k -> 1.0 / Math.pow(k, 2.0), 1.644934066848226436),
			new InfiniteSeries(k -> 1.0 / Math.pow(k, 3.0), 1.202056903159594),
			new InfiniteSeries(k -> (Math.sqrt(k) + 1.0) / Math.pow(k, 2), 4.257309415533714),
			new InfiniteSeries(k -> Math.pow(k + Math.exp(1.0 / k), -Constants.SQRT2), 1.7137967355403014865),
			new InfiniteSeries(k -> Math.log(1 + k) / Math.pow(1 + k, 1.5), 3.932239737431101),
			new InfiniteSeries(k -> Math.log(1 + k) / Math.pow(1 + k, 2.0), 0.937548254315843),
			new InfiniteSeries(k -> 1.0 / (1.0 + k) / Math.pow(Math.log(1.0 + k), 2), 2.10974280123689197448),
			new InfiniteSeries(k -> Math.log((k + 1.0) / (k + 0.0)) * Math.log((k + 2.0) / (k + 1.0)),
					0.68472478856315712) };

	public static final void main(String[] args) {
		InfiniteSeries[][] ALL_SERIES = { ALTERNATING_SERIES, IRREGULAR_SERIES, LINEAR_SERIES, LOGARITHMIC_SERIES };
		String[] disp = { "Alternating", "Irregular", "Linear", "Logarithmic" };
		int i = 0;
		for (InfiniteSeries[] series_of_type : ALL_SERIES) {
			System.out.println(disp[i]);
			for (InfiniteSeries series : series_of_type) {
				Function<Long, Double> transformed;
				final SeriesAlgorithm toalt;
				if (i == 3) {
					toalt = new Ensemble(1e-7, 100, 5, 0);
					transformed = toalt.alternatingSeries(series.mySeries, 1L);
				} else {
					toalt = null;
					transformed = series.mySeries;
				}
				final SeriesAlgorithm alg = new Ensemble(1e-6, 1000, 10);
				final double limit = alg.limit(transformed, 1L, true);
				final double actual = series.myLimit;
				final double error = Math.abs(limit - actual);
				System.out.println(error);
				final int sev;
				if (i == 3) {
					sev = toalt.countEvaluations();
				} else {
					sev = alg.countEvaluations();
				}
				System.out.println(sev);
			}
			i += 1;
		}
	}

	private SeriesTest() {
	}
}
