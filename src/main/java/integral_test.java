import java.util.function.Function;

import integral.Quadrature;
import integral.Romberg.RombergExtrapolationMethod;
import integral.ClenshawCurtis.ClenshawCurtisExtrapolationMethod;
import integral.*;

public class integral_test {

	public static void main(final String[] args) {
		Function<Double, Double> func = x -> Math.sin(x);
		double a = 0.;
		double b = 1.;
		double eps = 1e-6;

		Quadrature[] quads = { new BulirschStoer(eps), new ClenshawCurtis(eps, ClenshawCurtisExtrapolationMethod.HAVIE),
				new ClenshawCurtis(eps, ClenshawCurtisExtrapolationMethod.OLIVER), new GaussKronrod(eps),
				new GaussLegendre(eps), new GaussLobatto(eps), new NewtonCotes(eps), new RmsRule(eps),
				new Romberg(eps, RombergExtrapolationMethod.RICHARDSON),
				new Romberg(eps, RombergExtrapolationMethod.CADRE), new Romberg(eps, RombergExtrapolationMethod.HAVIE),
				new Simpson(eps), new TanhSinh(eps), new Trapezoid(eps) };

		for (Quadrature quad : quads) {
			System.out.println(
					quad.getClass().getName() + "\t" + quad.integrate(func, a, b) + "\t" + quad.getEvaluations());
		}
	}
}
