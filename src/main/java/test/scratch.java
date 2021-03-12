package test;

import java.util.function.Function;
import integral.*;
import series.dp.*;
import utils.Sequences;

public class scratch {

    public static final void main(String[] args) {

	// define the integrand
	Function<Double, Double> f = x -> Math.cos(x * x);

	// split the integration region into subregions
	Iterable<Double> roots = Sequences.toIterable(k -> Math.sqrt(Math.PI * k), 0L);
	Quadrature integrator = new GaussKronrod(1e-8, 1000);
	Iterable<Double> integrals = integrator.integratePiecewise(f, roots);
	
	// evaluate the integral as an alternating series
	SeriesAlgorithm accelerator = new Cohen(1e-8, 1000, 1);
	System.out.println(accelerator.limit(integrals, true));
    }
}
