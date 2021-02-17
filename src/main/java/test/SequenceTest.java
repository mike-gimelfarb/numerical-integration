package test;

import java.util.Iterator;
import java.util.function.Function;

import series.dp.*;
import utils.Sequences;

public final class SequenceTest {

    public static final class Sequence {

	public final Iterable<Double> mySequence;
	public final double myLimit;

	public Sequence(final Iterable<Double> sequence, final double limit) {
	    mySequence = sequence;
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

    public static final Iterable<Double> firstOrderRecurrence(final Function<Double, Double> f, final double a0) {
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

    public static final Sequence[] EXPLICIT = {

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

    public static final Sequence[] RECURRENCE = {

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

    public static final void main(String[] args) {
	Sequence[][] ALL_SEQ = { EXPLICIT, RECURRENCE };
	String[] disp = { "Explicit", "Recurrence" };
	int i = 0;
	for (Sequence[] seq_of_type : ALL_SEQ) {
	    System.out.println("\n" + disp[i]);
	    for (Sequence seq : seq_of_type) {
		final double tol = 1e-8;
		final int nconv = 5;
		final Iterable<Double> func = seq.mySequence;
		final SeriesAlgorithm alg = new Ensemble(tol, 500, nconv, 0);
		final double limit = alg.limit(func, false);
		final double actual = seq.myLimit;
		final double error = Math.abs(limit - actual);
		final int sev = alg.countEvaluations();
		System.out.println(error);
		System.out.println(sev);
	    }
	    ++i;
	}
    }

    private SequenceTest() {
    }
}
