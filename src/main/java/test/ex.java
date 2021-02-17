package test;

import java.util.function.Function;

import series.dp.*;
import utils.Sequences;

public class ex {

    public static final void main(String[] args) {
	Function<Long, Double> f = k -> Math.pow(1 + 1. / k, k);
	Iterable<Double> seq = Sequences.toIterable(f, 1L);
	SeriesAlgorithm alg = new Ensemble(1e-6, 100, 2, 1);
	System.out.println(alg.limit(seq, false));
	System.out.println(alg.countEvaluations());
    }
}
