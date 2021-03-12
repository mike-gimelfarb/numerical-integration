package test;

import java.util.Iterator;
import java.util.function.Function;

import utils.Constants;
import utils.Sequences;

public final class TestProblems {

    // **************************************************************************
    // SEQUENCES
    // **************************************************************************
    public static final class Sequence {

	public final Iterable<Double> mySequence;
	public final double myLimit;

	public Sequence(final Iterable<Double> sequence, final double limit) {
	    mySequence = sequence;
	    myLimit = limit;
	}
    }

    private static final Iterable<Double> firstOrderRecurrence(final Function<Double, Double> f, final double a0) {
	return () -> new Iterator<>() {

	    double a = a0;

	    @Override
	    public boolean hasNext() {
		return true;
	    }

	    @Override
	    public Double next() {
		a = f.apply(a);
		return a;
	    }
	};
    }

    private static final double logFactorial(long n) {
	if (n <= 1) {
	    return 0.0;
	} else {
	    return Math.log(n) + logFactorial(n - 1);
	}
    }

    public static final Sequence[] EXPLICIT_SEQUENCE = {

	    // functions of the form a^g(k)
	    new Sequence(Sequences.toIterable(n -> Math.pow(0.5, 1.0 / n), 2), 1.0),
	    new Sequence(Sequences.toIterable(n -> Math.pow(2.0, 1.0 / n), 2), 1.0),

	    // functions of the form k^g(k)
	    new Sequence(Sequences.toIterable(n -> Math.exp(Math.log(n) / n), 2), 1.0),

	    // functions of the form f(k)^g(k)
	    new Sequence(Sequences.toIterable(n -> Math.pow(1.0 + 1.0 / n, n), 2), Math.E),
	    new Sequence(Sequences.toIterable(n -> Math.pow(1.0 - 1.0 / n, n), 2), 1.0 / Math.E),

	    // sums, products and composites
	    new Sequence(Sequences.toIterable(n -> (Math.pow(2, 1.0 / n) - 1.0) / (1.0 / n), 2), Math.log(2)),
	    new Sequence(Sequences.toIterable(n -> (Math.exp(1.0 / n) - 1.0) / (1.0 / n), 2), 1.0),

	    // logs
	    new Sequence(Sequences.toIterable(n -> Math.log(1.0 + 1.0 / n) / (1.0 / n), 2), 1.0),
	    new Sequence(Sequences.toIterable(n -> Math.log(n) / n, 2), 0.0),

	    // trig
	    new Sequence(Sequences.toIterable(n -> Math.sin(1.0 / n) / (1.0 / n), 2), 1.0),
	    new Sequence(Sequences.toIterable(n -> Math.sin(5.0 / n) / (1.0 / n), 2), 5.0),
	    new Sequence(Sequences.toIterable(n -> n * Math.sin(1.0 / n), 2), 1.0),
	    new Sequence(Sequences.toIterable(n -> (1.0 - Math.cos(1.0 / n)) / (1.0 / n), 2), 0.0),
	    new Sequence(Sequences.toIterable(n -> (1.0 - Math.cos(1.0 / n)) / (1.0 / n / n), 2), 0.5),

	    // special
	    new Sequence(Sequences.toIterable(n -> n / Math.exp(logFactorial(n) / n), 2), Math.E), };

    public static final Sequence[] RECURRENCE_SEQUENCE = {

	    // functional
	    new Sequence(firstOrderRecurrence(x -> (x + 2.) / (3 * x + 2.), 2.), 2. / 3),
	    new Sequence(firstOrderRecurrence(x -> Math.sqrt(2. * x - 1.), 2.), 1.),
	    new Sequence(firstOrderRecurrence(x -> 0.5 * (1. + x * x * x), 0.), 0.61803398875),
	    new Sequence(firstOrderRecurrence(x -> 0.5 * (1. + x * x), 0.), 1.),

	    // trig
	    new Sequence(firstOrderRecurrence(x -> Math.cos(x), 0.), 0.739085133215),
	    new Sequence(firstOrderRecurrence(x -> Math.sin(x + 1), 0.), 0.934563210752),

	    // logistic map
	    new Sequence(firstOrderRecurrence(x -> 0.5 * x * (1 - x), .5), 0.),
	    new Sequence(firstOrderRecurrence(x -> 1.5 * x * (1 - x), .5), 1. / 3),
	    new Sequence(firstOrderRecurrence(x -> 2.5 * x * (1 - x), .5), 0.6),

	    // lambert W
	    new Sequence(firstOrderRecurrence(x -> Math.log(4. / x), 0.5), 1.20216787319704),
	    new Sequence(firstOrderRecurrence(x -> Math.log(3. / x), 0.5), 1.04990889496404),
	    new Sequence(firstOrderRecurrence(x -> Math.exp(-x), 0.), 0.567143290409783),
	    new Sequence(firstOrderRecurrence(x -> Math.exp(-x + 0.5), 0.), 0.76624860816175),
	    new Sequence(firstOrderRecurrence(x -> Math.exp(-x + 0.75), 0.), 0.87898614436894),
	    new Sequence(firstOrderRecurrence(x -> Math.exp(-x + 1), 0.), 1.0) };

    // **************************************************************************
    // INFINITE SERIES
    // **************************************************************************
    public static final class InfiniteSeries {

	public final Function<Long, Double> mySeries;
	public final double myLimit;

	public InfiniteSeries(final Function<Long, Double> series, final double limit) {
	    mySeries = series;
	    myLimit = limit;
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

    public static final InfiniteSeries[] FAST_CONVERGING_SERIES = {

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

    public static final InfiniteSeries[] SLOW_CONVERGING_SERIES = {

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
    public static final InfiniteSeries[] ARCANE_SERIES = {

	    // alternating
	    new InfiniteSeries(k -> Math.pow(-1, k) * Math.abs(Math.sin(k)) / k, -0.4050635),
	    new InfiniteSeries(k -> Math.pow(-1, k) / (k + Math.sin(k)), -0.3545469),

	    // really weird
	    new InfiniteSeries(k -> 1.0 / (Math.pow(k, 3) * Math.sin(k * Math.PI * Constants.SQRT2)), -0.791727),
	    new InfiniteSeries(k -> Math.sin(1 + k) / Math.log(1 + k), 0.6839137) };

    // **************************************************************************
    // INTEGRATION
    // **************************************************************************
    public static final class Integral {

	public final Function<Double, Double> myF;
	public final double myLower, myUpper;
	public final double myValue;

	public Integral(final Function<Double, Double> f, final double a, final double b, final double value) {
	    myF = f;
	    myLower = a;
	    myUpper = b;
	    myValue = value;
	}
    }

    public static final Integral[] SIMPLE_INTEGRALS = {

	    // polynomials
	    new Integral(x -> x * x - 2 * x + 3, 0, 1, 2.33333333333333), //
	    new Integral(x -> x * x * x - 2 * x * x + 3 * x - 4, 0, 1, -2.91666666666666), //

	    // rationals
	    new Integral(x -> x <= 0 ? 0 : x / (Math.exp(x) - 1), 0, 1, 0.7775046341122), //
	    new Integral(x -> 1 / (1 + x), 0, 1, 0.6931471805599), //
	    new Integral(x -> 1 / (1 + x * x), 0, 1, 0.78539816339744), //
	    new Integral(x -> 1 / (1 + Math.pow(x, 4)), 0, 1, 0.866972987), //
	    new Integral(x -> 1 / (1 + Math.exp(x)), 0, 1, 0.3798854930417), //
	    new Integral(x -> 1 / (Math.pow(x, 4) + x * x + 0.9), -1, 1, 1.5822329637296729), //
	    new Integral(x -> 1 / (1 + Math.sqrt(Math.max(x, 0.))), 0, 1, 0.6137056388801), //
	    new Integral(x -> 1 / Math.pow(4 * x * x + 1, 2), -1, 1, 0.753574358897045), //
	    new Integral(x -> 1 / (1 - 0.5 * Math.pow(x, 4)), 0, 1, 1.14366725406941), //
	    new Integral(x -> 1 / (1 - 0.98 * Math.pow(x, 4)), 0, 1, 1.896335631177699), //
	    new Integral(x -> 1 / (1 - 0.998 * Math.pow(x, 4)), 0, 1, 2.4670706247423), //
	    new Integral(x -> 4 / (1 + 256 * Math.pow(x - 3. / 8, 2)), 0, 1, 0.719193830921), //

	    // logarithms
	    new Integral(x -> x <= 0 ? 0 : x * Math.log(x), 0, 1, -0.25), //
	    new Integral(x -> x <= 0 ? 0 : x * x * Math.log(x), 0, 1, -0.1111111111111), //
	    new Integral(x -> x <= 0 ? 0 : x * Math.pow(Math.log(x), 2), 0, 1, 0.25), //
	    new Integral(x -> x * Math.log(1 + x), 0, 1, 0.25), //
	    new Integral(x -> x <= 0 ? 0 : Math.log(1 + x) / x, 0, 1, 0.82246703342411), //
	    new Integral(x -> Math.pow(x, 4) * Math.log(x + Math.sqrt(x * x + 1)), 0, 2, 8.1533641198111), //

	    // others
	    new Integral(x -> x * x * Math.atan(x), 0, 1, 0.2106572512258), //
	    new Integral(x -> Math.exp(x) * Math.cos(x), 0, Math.PI / 2, 1.9052386904826), //
	    new Integral(x -> 23. / 25 * Math.cosh(x) - Math.cos(x), -1, 1, 0.4794282266888), //
	    new Integral(x -> Math.atan(Math.sqrt(x * x + 1)) / Math.pow(x * x + 1, 1.5), 0, 1, 0.590489270886385), //
	    new Integral(x -> Math.atan(Math.sqrt(x * x + 2)) / (1 + x * x) / Math.sqrt(2 + x * x), 0, 1, 0.51404189589) //
    };

    public static final Integral[] OSCILLATING_INTEGRALS = {

	    // rationals
	    new Integral(x -> 2. / (2 + Math.sin(10 * Math.PI * x)), 0, 1, 1.15470053837925), //

	    // oscillating
	    new Integral(x -> x < 0 ? 0 : Math.sqrt(10000 * Math.PI * Math.PI - x * x) * Math.sin(x), 0, 100 * Math.PI,
		    298.43571649436), //
	    new Integral(x -> Math.exp(-Math.sqrt(Math.max(x, 0))) * Math.sin(x), 0, 2 * Math.PI, 0.38749284555451), //
	    new Integral(x -> Math.cos(32 * Math.sin(x)), 0, Math.PI, 0.4337880026347), //
	    new Integral(x -> Math.cos(64 * Math.sin(x)), 0, Math.PI, 0.2908801021737), //
	    new Integral(x -> Math.cos(Math.cos(x) + 3 * Math.sin(x) + 2 * Math.cos(2 * x) + 3 * Math.cos(3 * x)), 0,
		    Math.PI, 0.29101878286), //

	    // rapidly oscillating
	    new Integral(x -> x * Math.pow(Math.cos(20 * x), 2), 0, Math.PI, 2.467401100272339), //
	    new Integral(x -> x * Math.sin(30 * x) * Math.cos(x), 0, 2 * Math.PI, -0.20967247966116), //
	    new Integral(x -> x * Math.sin(30 * x) * Math.cos(50 * x), 0, 2 * Math.PI, 0.117809724509617), //
	    new Integral(
		    x -> x >= 2 * Math.PI ? 0 : x * Math.sin(30 * x) / Math.sqrt(1. - x * x / (4. * Math.PI * Math.PI)),
		    0, 2 * Math.PI, -2.54325961889353), //
	    new Integral(x -> x == 0 ? 0 : 50 * Math.pow(Math.sin(50 * Math.PI * x) / (50 * Math.PI * x), 2), 0, 1,
		    0.49898680869304), //
	    new Integral(x -> Math.cos(1. / x) / x, 0.01, 1, -0.3425527480435786), //
    };

    public static final Integral[] IRREGULAR_INTEGRALS = { //

	    // vertical asymptotes at endpoints
	    new Integral(x -> Math.sqrt(Math.max(x, 0)), 0, 1, 0.66666666666666), //
	    new Integral(x -> Math.pow(Math.max(x, 0), 0.25), 0, 1, 0.8), //
	    new Integral(x -> Math.pow(Math.max(x, 0), 0.125), 0, 1, 0.88888888888888), //
	    new Integral(x -> x <= 0 ? 0 : Math.sqrt(x) * Math.log(x), 0, 1, -4. / 9), //
	    new Integral(x -> Math.sqrt(1 - x * x), -1, 1, Math.PI / 2), //
	    new Integral(x -> Math.sqrt(Math.max(x, 0)) * Math.abs(Math.cos(6 * Math.PI * x)), 0, 1, 0.423474113339959), //

	    // spikes
	    new Integral(x -> Math.sqrt(Math.abs(x + 0.5)), -1, 1, 1.4604471317871), //
	    new Integral(x -> 1 / (Math.exp(-Math.abs(2 * x)) + 1), -1, 1, 1.433780830483), //
	    new Integral(x -> 1 / (Math.exp(-Math.abs(2 * Math.sin(2 * Math.PI * x))) + 1), -1, 1, 1.5288524565498178), //
	    new Integral(x -> Math.log(Math.abs((Math.tan(x) + Math.sqrt(7)) / (Math.tan(x) - Math.sqrt(7)))),
		    Math.PI / 3, Math.PI / 2, 0.8889149278163532), //

	    // plateau and sharp edges or other fine features
	    new Integral(x -> Math.exp(x) * Math.sin(3 * x) * Math.tanh(5 * Math.cos(30 * x)), -1, 1, -0.017790593076), //
    };

    private TestProblems() {
    }
}
