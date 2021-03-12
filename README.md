# Summary
This library contains standard and state-of-the-art algorithms for estimating integrals of functions in 1D and limits of sequences and series in Java, many of which are not currently implemented in Commons Math.

The main features are:
- over a dozen algorithms for estimating integrals of 1D functions
- support for variable transformations, oscillatory and indefinite integrals
- algorithms for estimating limits of sequences and series in double and high (e.g. arbitrary) precision

# Examples

## Limits and Series -- Simple Example
The framework currently supports many algorithms for accelerating convergence of infinite series, including methods of Aitken, Cohen, Theta (Brezinski and iterated), Levin (D, T, U and V transforms), Richardson, Wynn's Epsilon and Wynn's Rho, and an ensemble method.

Consider the following infinite series:

<a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{k=1}^{\infty}{\frac{1}{k^2}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{k=1}^{\infty}{\frac{1}{k^2}}" title="\sum_{k=1}^{\infty}{\frac{1}{k^2}}" /></a>

First, import the necessary packages:

```java
import java.util.function.Function;
import series.dp.*;
import utils.Sequences;
```

The Ensemble algorithm is the best bet if additional information about the series is not provided:

```java
// define the series terms
Function<Long,Double> f = k -> 1.0 / Math.pow(k, 2);
Iterable<Double> iter = Sequences.toIterable(f, 1);

// find limit
SeriesAlgorithm alg = new Ensemble(1e-8, 1000, 5);
System.out.println(alg.limit(iter, true));
```

The estimate of the limit is 1.64493406, the exact value is 1.64493406 (to 8 decimal places).

## Limits and Series -- Complex Example
Sometimes existing methods fail if the series is very slowly converging. Consider the following very difficult series:

<a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{k=2}^{\infty}{\frac{1}{k&space;(\log{k})^2}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{k=2}^{\infty}{\frac{1}{k&space;(\log{k})^2}}" title="\sum_{k=2}^{\infty}{\frac{1}{k (\log{k})^2}}" /></a>

Applying the Ensemble algorithm directly will return 2.06035757 which is not the correct limit. However, using a trick of van Wijngaarden, we can convert the series to an alternating series with the same limit:

```java
// define the series terms
Function<Long,Double> f = k -> 1.0 / k / Math.pow(Math.log(k), 2);
	
// create transform of original series
Function<Long, Double> f2 = (new Ensemble(1e-8, 1000, 5)).toAlternatingSeries(f, 2);
Iterable<Double> iter = Sequences.toIterable(f2, 1);
	
// find limit of transformed series
SeriesAlgorithm alg = new Ensemble(1e-8, 1000, 5);
System.out.println(alg.limit(iter, true));
```

Now, the estimate of the limit is 2.10972103.

## Definite Integrals -- Simple Example
The framework is also very capable of numerically evaluating integrals of 1d functions, including adaptive Clenshaw-Curtis, adaptive Gaussian (Gauss-Kronrod, Gauss-Legendre, Gauss-Lobatto, RMS improvement of Gauss-Kronrod), Newton-Cotes, Romberg, adaptive Simpson (global and local versions), and Tanh-Sinh quadrature schemes.

Consider the following definite integral:

<a href="https://www.codecogs.com/eqnedit.php?latex=\int_{-100}^{100}\frac{1}{1&plus;x^2}\,\mathrm{d}x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\int_{-100}^{100}\frac{1}{1&plus;x^2}\,\mathrm{d}x" title="\int_{-100}^{100}\frac{1}{1+x^2}\,\mathrm{d}x" /></a>

First, import the necessary packages:

```java
import java.util.function.Function;
import integral.*;
```

We will use Gaussian quadrature (Gauss-Kronrod) to evaluate this integral:

```java
// define the integrand
Function<Double, Double> f = x -> 1.0 / (1.0 + x * x);
	
// integrate
Quadrature integrator = new GaussKronrod(1e-8, 1000);
System.out.println(integrator.integrate(f, -100.0, 100.0));
```

The predicted value is 3.12159332, which is correct to 8 decimal places.

## Definite Integrals -- Complex Example

The following integral is improper, and is also highly oscillatory:

<a href="https://www.codecogs.com/eqnedit.php?latex=\int_{0}^{\infty}\cos(x^2)\,\mathrm{d}x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\int_{0}^{\infty}\cos(x^2)\,\mathrm{d}x" title="\int_{0}^{\infty}\cos(x^2)\,\mathrm{d}x" /></a>

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
