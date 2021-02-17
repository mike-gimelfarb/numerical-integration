package test;

import java.util.function.Function;

import series.dp.*;
import utils.Constants;
import utils.Sequences;

public final class SeriesTest {

    public static final class InfiniteSeries {

	public final Function<Long, Double> mySeries;
	public final double myLimit;

	public InfiniteSeries(final Function<Long, Double> series, final double limit) {
	    mySeries = series;
	    myLimit = limit;
	}
    }

    public static final double logFactorial(long n) {
	if (n <= 1) {
	    return 0.0;
	} else {
	    return Math.log(n) + logFactorial(n - 1);
	}
    }

    public static final InfiniteSeries[] ALTERNATING_SERIES = {

	    // reciprocal power
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / Math.pow(k, 1.5), 0.76514702462540794),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / k, Math.log(2)),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / (2 * k - 1.), Math.PI / 4),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / Math.pow(k, 0.5), 0.60489864342163037),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / Math.pow(k, 0.25), 0.5544873859140731),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / Math.pow(k, 1.0 + 1.0 / k), 0.77951153739328),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) * Math.sin(1.0 / k), 0.550796848133929),

	    // slow converging
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / Math.log(k + 1), 0.9242998972229),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) * Math.log(k + 1) / (k + 1), 0.1598689037),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / Math.sqrt(Math.log(k + 1)), 0.69024445098278),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) / Math.log(Math.log(k + 1)), -11.4779680139871),
	    new InfiniteSeries(k -> Math.pow(-1, k - 1) * Math.log(k + 1) / Math.sqrt(k + 1), 0.1932888316) };

    public static final InfiniteSeries[] IRREGULAR_SERIES = {
	    new InfiniteSeries(k -> Math.sin(k) / k / Math.pow(2, k), 0.522937836179332),
	    new InfiniteSeries(k -> Math.sin(5 * k) / k / Math.pow(2, k), -0.509500940834413),
	    new InfiniteSeries(k -> Math.sin(k) / k, 1.0707963267948966),
	    new InfiniteSeries(k -> Math.sin(5 * k) / k, -0.92920367320510338),
	    new InfiniteSeries(k -> Math.pow(Math.sin(k), 3) / k, 0.7853981634080),
	    new InfiniteSeries(k -> Math.pow(-1, k / 2) / k, 0.43882457311747),
	    new InfiniteSeries(k -> Math.cos(1. * k) / k, 0.042019505825),
	    new InfiniteSeries(k -> Math.cos(5. * k) / k, -0.179771882642) };

    public static final InfiniteSeries[] FAST_CONVERGING = {

	    // geometric series
	    new InfiniteSeries(k -> Math.pow(0.8, k) / k, 1.609437912434100),
	    new InfiniteSeries(k -> Math.pow(0.995, k - 1), 200.),
	    new InfiniteSeries(k -> Math.pow(0.4, k - 1) + Math.pow(0.8, k - 1), 6.666666666666666),
	    new InfiniteSeries(k -> Math.pow(0.8, k - 1) * Math.pow(k - 1, 2), 180),

	    // exponential function
	    new InfiniteSeries(k -> 1.0 / Math.exp(logFactorial(k - 1)), Math.E),
	    new InfiniteSeries(k -> Math.exp((k - 1) * Math.log(2.) - logFactorial(k - 1)), Math.exp(2.)),

	    // others
	    new InfiniteSeries(k -> Math.sin(1.0 / k) / k, 1.472828231956185296295),
	    new InfiniteSeries(k -> Math.sin(1.0 / k) / Math.sqrt(k), 2.432932908227431),
	    new InfiniteSeries(k -> Math.exp(2 * logFactorial(k - 1) - logFactorial(2 * (k - 1))), 1.736399858718715) };

    public static final InfiniteSeries[] SLOW_CONVERGING = {

	    // zeta function
	    new InfiniteSeries(k -> 1 / Math.pow(k, 2), 1.644934066848226436),
	    new InfiniteSeries(k -> 1 / Math.pow(k, 1.5), 2.612375348685488343),
	    new InfiniteSeries(k -> 1 / Math.pow(k, 1.25), 4.59511182584294338),
	    new InfiniteSeries(k -> 1 / Math.pow(k, 1.125), 8.5862412945105753),

	    // rationals and powers
	    new InfiniteSeries(k -> (2 * k - 1.) / (k * (k + 1.) * (k + 2.)), 0.75),
	    new InfiniteSeries(k -> (Math.sqrt(k) + 1) / Math.pow(k, 2), 4.257309415533714),
	    new InfiniteSeries(k -> 1.0 / Math.pow(2, Math.sqrt(k)), 3.788219230647953951),
	    new InfiniteSeries(k -> Math.pow(k + Math.exp(1. / k), -Constants.SQRT2), 1.713796735540301),
	    new InfiniteSeries(k -> Math.pow(k, 1.0 / k / k) - 1., 0.97149903428330),

	    // logs
	    new InfiniteSeries(k -> k == 1 ? 1. : 1. / k + Math.log(1 - 1. / k), 0.577215664901532),
	    new InfiniteSeries(k -> Math.log((k + 1.) / k) * Math.log((k + 2.) / (k + 1.)), 0.68472478856315712),
	    new InfiniteSeries(k -> Math.log(k) / Math.pow(k, 2), 0.937548254315843),
	    new InfiniteSeries(k -> Math.log(k) / Math.pow(k, 1.5), 3.932239737431101),
	    new InfiniteSeries(k -> Math.log(k) / Math.pow(k, 1.25), 15.92966500287924),
	    new InfiniteSeries(k -> Math.log(k) / Math.pow(k, 2), 0.937548254315843),
	    new InfiniteSeries(k -> Math.pow(Math.log(k) / k, 2), 1.9892802342989010),
	    new InfiniteSeries(k -> Math.pow(Math.log(k), 2) / Math.pow(k, 1.25), 127.989866745379),
	    new InfiniteSeries(k -> 1. / ((k + 1.) * Math.pow(Math.log(k + 1), 2)), 2.109742801236892),
	    new InfiniteSeries(k -> 1.0 / Math.pow(Math.log(k + 1), Math.log(k + 1)), 5.716970612990),
	    new InfiniteSeries(k -> 1.0 / Math.pow(3, Math.log(k)), 10.7250727422677186),
	    new InfiniteSeries(k -> Math.log(1. + 1.0 / k) / k, 1.25774688694436963),
	    new InfiniteSeries(k -> Math.log(1. + 1.0 / k) / Math.sqrt(k), 2.166970206923584479),

	    // misc
	    new InfiniteSeries(k -> Math.sin(1. / k) / Math.sqrt(k), 2.4329329082274314),
	    new InfiniteSeries(k -> Math.sin(1. / k) * Math.log(Math.cos(1 / Math.sqrt(k))), -0.852090754198727956015),
	    new InfiniteSeries(k -> Math.log(k * Math.sin(1. / k)), -0.280556336229) };

    // these series are too much for some solvers
    public static final InfiniteSeries[] ARCANE = {

	    // alternating
	    new InfiniteSeries(k -> Math.pow(-1, k) * Math.abs(Math.sin(k)) / k, -0.4050635),
	    new InfiniteSeries(k -> Math.pow(-1, k) / (k + Math.sin(k)), -0.3545469),

	    // really weird
	    new InfiniteSeries(k -> 1.0 / (Math.pow(k, 3) * Math.sin(k * Math.PI * Constants.SQRT2)), -0.791727),
	    new InfiniteSeries(k -> Math.sin(1 + k) / Math.log(1 + k), 0.6839137) };

    public static final void main(String[] args) {
	InfiniteSeries[][] ALL_SERIES = { ALTERNATING_SERIES, IRREGULAR_SERIES, FAST_CONVERGING, SLOW_CONVERGING,
		ARCANE };
	String[] disp = { "Alternating", "Irregular", "Fast", "Slow", "Arcane" };
	int i = 0;
	for (InfiniteSeries[] series_of_type : ALL_SERIES) {
	    System.out.println("\n" + disp[i]);
	    for (InfiniteSeries series : series_of_type) {

		final double tol;
		final int nconv;
		if (i == 1 || i == 4) {
		    tol = 1e-6;
		    nconv = 20;
		} else {
		    tol = 1e-8;
		    nconv = 5;
		}
		final Function<? super Long, Double> func;
		final SeriesAlgorithm algtoalt;
		if (i == 3) {
		    algtoalt = new Ensemble(tol, 50, 2, 0);
		    func = algtoalt.toAlternatingSeries(series.mySeries, 1L);
		} else {
		    algtoalt = null;
		    func = series.mySeries;
		}
		final SeriesAlgorithm alg = new Ensemble(tol, 500, nconv, 0);
		final double limit = alg.limit(Sequences.toIterable(func, 1L), true, 3);
		final double actual = series.myLimit;
		final double error = Math.abs(limit - actual);
		System.out.println(error);
		final int sev;
		if (algtoalt != null) {
		    sev = algtoalt.countEvaluations();
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
