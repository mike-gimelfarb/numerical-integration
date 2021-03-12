package test;

import java.text.DecimalFormat;

import integral.*;
import integral.ClenshawCurtis.ClenshawCurtisExtrapolationMethod;
import integral.Quadrature.QuadratureResult;
import integral.Romberg.RombergExtrapolationMethod;
import integral.Simpson.SimpsonAdaptationType;
import test.TestProblems.Integral;

public final class FiniteIntegralTest {

    private static final String pad(final String str, final int len) {
	if (str.length() >= len) {
	    return str;
	}
	String result = str;
	for (int i = str.length() + 1; i <= len; ++i) {
	    result += " ";
	}
	return result;
    }

    public static final void main(String[] args) {

	// algorithms
	final double tol = 1e-8;
	final int maxev = 30000;
	Quadrature[] ALGOS = { //
		new ClenshawCurtis(tol, maxev, ClenshawCurtisExtrapolationMethod.HAVIE), //
		new ClenshawCurtis(tol, maxev, ClenshawCurtisExtrapolationMethod.OLIVER), //
		new GaussKronrod(tol, maxev), //
		new GaussLegendre(tol, maxev), //
		new GaussLobatto(tol, maxev), //
		new NewtonCotes(tol, maxev), //
		new RmsRule(tol, maxev), //
		new Romberg(tol, maxev, RombergExtrapolationMethod.CADRE), //
		new Romberg(tol, maxev, RombergExtrapolationMethod.HAVIE), //
		new Romberg(tol, maxev, RombergExtrapolationMethod.RICHARDSON), //
		new Simpson(tol, maxev, SimpsonAdaptationType.GLOBAL), //
		new Simpson(tol, maxev, SimpsonAdaptationType.LOCAL), //
		new TanhSinh(tol, maxev) };

	Integral[][] ALL = { //
		TestProblems.SIMPLE_INTEGRALS, //
		TestProblems.OSCILLATING_INTEGRALS, //
		TestProblems.IRREGULAR_INTEGRALS };
	String[] disp = { "Simple", "Oscillating", "Irregular" };

	int nintegs = 0;
	for (final Integral[] integrals_of_type : ALL) {
	    nintegs += integrals_of_type.length;
	}
	final double[][] errors = new double[nintegs][ALGOS.length];
	final int[][] fev = new int[nintegs][ALGOS.length];
	int i = 0;
	int row = 0;
	for (Integral[] integral_of_type : ALL) {
	    for (Integral integral : integral_of_type) {
		int col = 0;
		for (Quadrature alg : ALGOS) {
		    final QuadratureResult result = alg.integrate(integral.myF, integral.myLower, integral.myUpper);
		    final double value = result.estimate;
		    final double actual = integral.myValue;
		    errors[row][col] = Math.abs(value - actual);
		    fev[row][col] = result.evaluations;
		    ++col;
		}
		++row;
	    }
	    ++i;
	}

	DecimalFormat format = new DecimalFormat("0.###E0");
	format.setMaximumFractionDigits(4);
	row = 0;
	i = 0;
	for (Integral[] integral_of_type : ALL) {
	    System.out.println("\n" + disp[i]);
	    for (final Quadrature alg : ALGOS) {
		System.out.print(pad(alg.getName(), 25));
	    }
	    System.out.println("");

	    for (Integral integral : integral_of_type) {
		for (int col = 0; col < ALGOS.length; ++col) {
		    String str = format.format(errors[row][col]) + " [" + fev[row][col] + "]";
		    System.out.print(pad(String.format("%20s", str), 25));
		}
		System.out.println("");
		++row;
	    }
	    i += 1;
	}
    }

    private FiniteIntegralTest() {
    }
}
