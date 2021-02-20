package test;

import java.util.function.Function;

import series.dp.*;
import test.TestProblems.InfiniteSeries;
import utils.Sequences;

public final class SeriesTest {

    public static final void main(String[] args) {
	InfiniteSeries[][] ALL_SERIES = { TestProblems.ALTERNATING_SERIES, TestProblems.IRREGULAR_SERIES,
		TestProblems.FAST_CONVERGING_SERIES, TestProblems.SLOW_CONVERGING_SERIES, TestProblems.ARCANE_SERIES };
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
