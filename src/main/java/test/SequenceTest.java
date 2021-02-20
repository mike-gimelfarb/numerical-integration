package test;

import series.dp.*;
import test.TestProblems.Sequence;

public final class SequenceTest {

    public static final void main(String[] args) {
	Sequence[][] ALL_SEQ = { TestProblems.EXPLICIT_SEQUENCE, TestProblems.RECURRENCE_SEQUENCE };
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
