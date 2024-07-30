# Summary
This library contains standard and state-of-the-art algorithms for estimating integrals of functions in 1D and limits of sequences and series in Java, many of which are not currently implemented in Commons Math.

The main features are:
- over a dozen algorithms for estimating integrals of 1D functions
- support for variable transformations, oscillatory and indefinite integrals
- algorithms for estimating limits of sequences and series in double and high (e.g. arbitrary) precision

# Examples

The necessary packages from the library to run the examples can be imported as follows

```java
import java.util.function.Function;
import integral.*;
import series.dp.*;
import utils.Sequences;
```

## Limits and Series -- Simple Example
The framework currently supports many algorithms for accelerating convergence of infinite series, including methods of Aitken, Cohen, Theta (Brezinski and iterated), Levin (D, T, U and V transforms), Richardson, Wynn's Epsilon and Wynn's Rho, and an ensemble method.

Consider the following infinite series:

$$\sum_{k=1}^{\infty}{\frac{1}{k^2}}$$

The Ensemble algorithm is the best bet if additional information about the series is not provided:

```java
// define the series terms
Iterable<Double> iter = Sequences.toIterable(k -> 1.0 / Math.pow(k, 2), 1);

// find limit
SeriesAlgorithm alg = new Ensemble(1e-8, 1000, 5);
System.out.println(alg.limit(iter, true));
```

The estimate of the limit is 1.64493406, the exact value is 1.64493406 (to 8 decimal places).

## Limits and Series -- Complex Example
Sometimes existing methods fail if the series is very slowly converging. Consider the following very difficult series:

$$\sum_{k=2}^{\infty}{\frac{1}{k (\log{k})^2}}$$

Applying the Ensemble algorithm directly will return 2.06035757 which is not the correct limit. However, using a trick of van Wijngaarden, we can convert the series to an alternating series with the same limit:

```java
// define the series terms
Function<Long, Double> f = k -> 1.0 / k / Math.pow(Math.log(k), 2);
	
// create transform of original series
Function<Long, Double> f2 = (new Ensemble(1e-8, 100, 2)).toAlternatingSeries(f, 2);
Iterable<Double> iter = Sequences.toIterable(f2, 1);
	
// find limit of transformed series
SeriesAlgorithm alg = new Ensemble(1e-8, 1000, 5);
System.out.println(alg.limit(iter, true));
```

Now, the estimate of the limit is 2.10973574, which is correct to 5 decimal places.

However, for terms that are sufficiently smooth and have an analytic continuation to the real numbers, it is possible to use the Euler-Maclaurin method to approximate the series by an imporoper integral, which often provides a more precise result. To run this, a method for evaluating integrals must be defined, how to break the improper integral up into proper integrals, and how to accelerate the resulting sequence of integral estimates.

```java
import series.special.EulerMaclaurin;

// define the series terms
Function<Double, Double> f = k -> 1.0 / k / Math.pow(Math.log(k), 2);
	
// find limit of transformed series
EulerMaclaurin algorithm = new EulerMaclaurin(1e-8, new GaussKronrod(1e-8, 1000), new Ensemble(1e-8, 1000, 5), k -> Math.pow(2, k));
System.out.println(algorithm.limit(f, 2L));
```

This time, the estimate of the limit is 2.10974280 and is correct to 8 decimals.

## Definite Integrals -- Simple Example
The framework is also very capable of numerically evaluating integrals of 1d functions, including adaptive Clenshaw-Curtis, adaptive Gaussian (Gauss-Kronrod, Gauss-Legendre, Gauss-Lobatto, RMS improvement of Gauss-Kronrod), Newton-Cotes, Romberg, adaptive Simpson (global and local versions), and Tanh-Sinh quadrature schemes.

Consider the following definite integral:

$$\int_{-100}^{100}\frac{1}{1&plus;x^2} \mathrm{d}x$$

We will use Gaussian quadrature (Gauss-Kronrod) to evaluate this integral:

```java
// define the integrand
Function<Double, Double> f = x -> 1.0 / (1.0 + x * x);
	
// integrate
Quadrature integrator = new GaussKronrod(1e-8, 1000);
System.out.println(integrator.integrate(f, -100.0, 100.0));
```

The predicted value is 3.12159332, which is correct to 8 decimal places.

## Definite Integrals -- Infinite Bounds Example

The following integral is improper and hard to evaluate directly:

$$\int_{-\infty}^{-1}\frac{e^x}{x} \mathrm{d}x$$

However, the package provides a transformation of variables to deal with integrals of this form:

$$[a,\infty]&space;:&space;t&space;\in&space;[0,1]&space;\to&space;(a&space;-&space;1)&space;&plus;&space;\frac{1}{t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?[a,\infty]&space;:&space;t&space;\in&space;[0,1]&space;\to&space;(a&space;-&space;1)&space;&plus;&space;\frac{1}{t}" title="[a,\infty] : t \in [0,1] \to (a - 1) + \frac{1}{t}" /></a>
<p></p>
<a href="https://www.codecogs.com/eqnedit.php?latex=[-\infty,b]&space;:&space;t&space;\in&space;[0,1]&space;\to&space;(b&space;&plus;&space;1)&space;-&space;\frac{1}{t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?[-\infty,b]&space;:&space;t&space;\in&space;[0,1]&space;\to&space;(b&space;&plus;&space;1)&space;-&space;\frac{1}{t}" title="[-\infty,b] : t \in [0,1] \to (b + 1) - \frac{1}{t}" /></a>

Integrals over [-infinity, infinity] are broken into two intervals at zero.

```java
// define the integrand
Function<Double, Double> f = x -> Math.exp(x) / x;

// split the integration region into subregions
Quadrature integrator = new TanhSinh(1e-8, 1000);
System.out.println(integrator.integrate(f, Double.NEGATIVE_INFINITY, -1.0));
```

The estimated value of the integral is -0.21938393, correct to 8 decimals.

## Definite Integrals -- Oscillatory with Infinite Bounds Example

The following integral is improper, and is also highly oscillatory:

$$\int_{0}^{\infty}\cos(x^2) \mathrm{d}x""

It is impossible to evaluate using any single integrator. However, it is possible to break up this integral over integration regions separated by the roots of the function, and estimate the above integral by evaluating an infinite series. The package provides a convenience method to do this, but requires the roots of the integrand to be fed to the method. The following example illustrates this:

```java
// define the integrand
Function<Double, Double> f = x -> Math.cos(x * x);

// split the integration region into subregions
Iterable<Double> roots = Sequences.toIterable(k -> Math.sqrt(Math.PI * k), 0L);
Quadrature integrator = new GaussKronrod(1e-8, 1000);
Iterable<Double> integrals = integrator.integratePiecewise(f, roots);
	
// evaluate the integral as an alternating series
SeriesAlgorithm accelerator = new Cohen(1e-8, 1000, 1);
System.out.println(accelerator.limit(integrals, true));
```

The estimated value of the integral is 0.62665707, which is correct to 8 decimal places.
