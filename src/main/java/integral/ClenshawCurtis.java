package integral;

import java.util.function.Function;

import utils.Constants;
import utils.SimpleMath;

/**
 * An adaptive numerical integrator based on the Clenshaw-Curtis quadrature
 * rule. There are two implementations:
 * <ol>
 * <li>Havie : the integral is approximated by Chebychev polynomials over each
 * subinterval, as introduced in [1]. This code is a translation of the Fortran
 * subroutine by John Burkardt.</li>
 * <li>Oliver : a doubly-adaptive Clenshaw-Curtis algorithm also using Chebychev
 * polynomials as introduced in [2].</li>
 * </ol>
 * 
 * <p>
 * References:
 * <ul>
 * <li>[1] Håvie, T. "On a modification of the Clenshaw-Curtis quadrature
 * formula." BIT Numerical Mathematics 9.4 (1969): 338-350.</li>
 * <li>[2] Oliver, J. "A doubly-adaptive Clenshaw-Curtis quadrature method." The
 * Computer Journal 15.2 (1972): 141-147.</li>
 * </ul>
 * </p>
 */
public final class ClenshawCurtis extends Quadrature {

	/**
	 * The extrapolation method to use for Clenshaw-Curtis integration.
	 */
	public static enum ClenshawCurtisExtrapolationMethod {
		HAVIE, OLIVER
	}

	private static final double[][] SIGMA = { { 0.455, 0.272, 0.606, 0.811, 0.908 },
			{ 0.550, 0.144, 0.257, 0.376, 0.511 }, { 0.667, 0.243, 0.366, 0.449, 0.522 },
			{ 0.780, 0.283, 0.468, 0.565, 0.624 }, { 0.855, 0.290, 0.494, 0.634, 0.714 },
			{ -1.00, 0.292, 0.499, 0.644, 0.745 } };

	private final int myMaxLevels;
	private final ClenshawCurtisExtrapolationMethod myMethod;

	/**
	 * Creates a new instance of the Clenshaw-Curtis integrator.
	 * 
	 * @param tolerance the smallest acceptable absolute change in integral
	 *                  estimates in consecutive iterations that indicates the
	 *                  algorithm has converged
	 * @param maxLevels determines the number of interpolation points
	 * @param method    the extrapolation method to use
	 */
	public ClenshawCurtis(final double tolerance, final int maxLevels, final ClenshawCurtisExtrapolationMethod method) {
		super(tolerance);
		myMaxLevels = maxLevels;
		myMethod = method;
	}

	/**
	 * Creates a new instance of the Clenshaw-Curtis integrator.
	 * 
	 * @param tolerance the smallest acceptable absolute change in integral
	 *                  estimates in consecutive iterations that indicates the
	 *                  algorithm has converged
	 * @param method    the extrapolation method to use
	 */
	public ClenshawCurtis(final double tolerance, final ClenshawCurtisExtrapolationMethod method) {
		this(tolerance, 12, method);
	}

	public ClenshawCurtis(final double tolerance) {
		this(tolerance, 12, ClenshawCurtisExtrapolationMethod.HAVIE);
	}

	@Override
	final double properIntegral(final Function<? super Double, Double> f, final double a, final double b) {
		switch (myMethod) {
		case HAVIE:
			return havie(f, a, b);
		case OLIVER:
			return oliver(f, a, b);
		default:
			return Double.NaN;
		}
	}

	private double havie(final Function<? super Double, Double> func, final double a, final double b) {
		final int nupper = myMaxLevels;
		final int len = 1 << nupper;
		final int hlen = len >> 1;
		final double[] acof = new double[hlen + 1], bcof = new double[hlen + 1], ccof = new double[len + 1];
		double a0, a1, a2, alf, alfnj, alfno, bet, betnj, betno, bounds, cof, cofmax, const1, const2, deln, deltan,
				error, etank, gamman, hnstep, one = 1.0, r1, r2, rk, rn, rnderr, rounde, tend, tnew, triarg, umid,
				wmean, xplus, xmin, xsink, epsout;
		int i, index, j, k, ksign, n, ncof, nhalf, nn, fev;

		// null integral
		if (a == b) {
			return 0.0;
		}

		// ROUNDE = RNDERR*(R1+R2*N), where R1, R2 are two empirical constants.
		// Set coefficients in formula for accumulated roundoff error.
		// N is the current number of function values used.
		rnderr = Constants.EPSILON;
		r1 = 1.0;
		r2 = 2.0;
		error = myTol;

		// Integration interval parameters.
		alf = 0.5 * (b - a);
		bet = 0.5 * (b + a);

		// Parameters for trigonometric recurrence relations.
		triarg = Math.atan(1.0);
		alfno = -1.0;

		// Parameters for integration stepsize and loops.
		rn = 2.0;
		n = 2;
		nhalf = 1;
		hnstep = 1.0;

		// Initial calculation for the end-point approximation.
		const1 = 0.5 * (func.apply(a) + func.apply(b));
		const2 = func.apply(bet);
		fev = 3;
		acof[1 - 1] = 0.5 * (const1 + const2);
		acof[2 - 1] = 0.5 * (const1 - const2);
		bcof[2 - 1] = acof[2 - 1];
		tend = 2.0 * (acof[1 - 1] - acof[2 - 1] / 3.0);

		// Start actual calculations.
		for (i = 1; i <= nupper; ++i) {

			// Compute function values.
			const1 = -Math.sin(triarg);
			const2 = 0.5 * alfno / const1;
			alfno = const1;
			betno = const2;
			gamman = 1.0 - 2.0 * alfno * alfno;
			deltan = -2.0 * alfno * betno;
			bcof[1 - 1] = 0.0;
			for (j = 1; j <= nhalf; ++j) {
				alfnj = gamman * const1 + deltan * const2;
				betnj = gamman * const2 - deltan * const1;
				xplus = alf * alfnj + bet;
				xmin = -alf * alfnj + bet;
				ccof[j - 1] = func.apply(xplus) + func.apply(xmin);
				fev += 2;
				bcof[1 - 1] += ccof[j - 1];
				const1 = alfnj;
				const2 = betnj;
			}
			bcof[1 - 1] *= 0.5 * hnstep;

			// Calculation of first B-coefficient finished compute the higher
			// coefficients if NHALF greater than one.
			if (nhalf > 1) {
				const1 = one;
				const2 = 0.0;
				ncof = nhalf - 1;
				ksign = -1;
				for (k = 1; k <= ncof; ++k) {

					// Compute trigonometric sum for B-coefficient.
					etank = gamman * const1 - deltan * const2;
					xsink = gamman * const2 + deltan * const1;
					cof = 2.0 * (2.0 * etank * etank - 1.0);
					a1 = 0.0;
					a0 = ccof[nhalf - 1];
					for (j = 1; j <= ncof; ++j) {
						a2 = a1;
						a1 = a0;
						index = nhalf - j;
						a0 = ccof[index - 1] + cof * a1 - a2;
					}
					bcof[k + 1 - 1] = hnstep * (a0 - a1) * etank;
					bcof[k + 1 - 1] = ksign * bcof[k + 1 - 1];
					ksign = -ksign;
					const1 = etank;
					const2 = xsink;
				}
			}

			// Calculation of B-coefficients finished.
			// Compute new modified mid-point approximation when the interval
			// of integration is divided in N equal sub intervals.
			umid = 0.0;
			rk = rn;
			nn = nhalf + 1;
			for (k = 1; k <= nn; ++k) {
				index = nn + 1 - k;
				umid += bcof[index - 1] / (rk * rk - one);
				rk -= 2.0;
			}
			umid *= -2.0;

			// Compute new C-coefficients for end-point approximation and largest
			// absolute value of coefficients.
			nn = n + 2;
			cofmax = 0.0;
			for (j = 1; j <= nhalf; ++j) {
				index = nn - j;
				ccof[j - 1] = 0.5 * (acof[j - 1] + bcof[j - 1]);
				ccof[index - 1] = 0.5 * (acof[j - 1] - bcof[j - 1]);
				const1 = Math.abs(ccof[j - 1]);
				cofmax = Math.max(cofmax, const1);
				const1 = Math.abs(ccof[index - 1]);
				cofmax = Math.max(cofmax, const1);
			}
			ccof[nhalf + 1 - 1] = acof[nhalf + 1 - 1];

			// Calculation of new coefficients finished.
			// Compute new end-point approximation when the interval of
			// integration is divided in 2N equal sub intervals.
			wmean = 0.5 * (tend + umid);
			bounds = 0.5 * (tend - umid);
			deln = 0.0;
			rk = 2.0 * rn;
			for (j = 1; j <= nhalf; ++j) {
				index = n + 2 - j;
				deln += ccof[index - 1] / (rk * rk - one);
				rk -= 2.0;
			}
			deln *= -2.0;
			tnew = wmean + deln;
			epsout = Math.abs(bounds / tnew);
			if (cofmax >= rnderr) {
				rounde = rnderr * (r1 + r2 * rn);
				epsout = Math.max(epsout, rounde);
				error = Math.max(error, rounde);
				if (error >= epsout) {

					// Required accuracy obtained or the maximum number of
					// function values used without obtaining the required accuracy.
					deln *= alf;
					myFEvals += fev;
					if (epsout <= myTol) {
						return alf * tnew;
					} else {
						return Double.NaN;
					}
				}
			}

			// If I = NUPPER then the required accuracy is not obtained.
			if (i == nupper) {

				// Required accuracy obtained or the maximum number of function
				// values used without obtaining the required accuracy.
				deln *= alf;
				myFEvals += fev;
				if (epsout <= myTol) {
					return alf * tnew;
				} else {
					return Double.NaN;
				}
			}
			System.arraycopy(ccof, 0, acof, 0, n);
			acof[n + 1 - 1] = ccof[n + 1 - 1];
			bcof[n + 1 - 1] = ccof[n + 1 - 1];
			tend = tnew;
			nhalf = n;
			n <<= 1;
			rn *= 2.0;
			hnstep *= 0.5;
			triarg *= 0.5;
		}
		myFEvals += fev;
		return Double.NaN;
	}

	private double oliver(final Function<? super Double, Double> f, double a, double b) {
		final double eta = 5.4210109e-20;
		final int divmax = 1 << myMaxLevels;
		double eps = myTol;
		int i, j = 0, m, mmax, n, n2, nmax, maxrule, order = 0, div, fev = 0;
		double c, cprev = 0.0, e = 0.0, eprev = 0.0, fmax = 0.0, fmin = 0.0, h = 0.0, hmin, iint = 0.0, intprev = 0.0,
				k, k1 = 0.0, re, x, xa, xb, xc, ans, error, acc = 0.;
		final double[] ec = new double[6], xs = new double[divmax], fs = new double[divmax], w1 = new double[129],
				fx = new double[129], t = new double[65], w = new double[126];
		boolean caution = false;

		ans = error = 0.0;
		if (a == b) {
			return ans;
		} else if (a > b) {
			c = b;
			b = a;
			a = c;
		}
		hmin = (b - a) / SimpleMath.pow(2.0, divmax);
		if (acc < eta) {
			acc = eta;
			acc *= 16.0;
		}
		x = 4;
		xa = 64;
		for (i = 1; i <= 6; ++i) {
			ec[i - 1] = xa / ((x * x - 1.0) * (x * x - 9.0));
			x *= 2.0;
			xa *= 2.0;
		}
		n = 4;
		n2 = 2;
		nmax = 128;
		m = mmax = nmax / n;
		t[0] = 1.0;
		t[nmax >> 1] = 0.0;
		maxrule = div = 0;
		w1[0] = -1.0;
		maxrule = quadrule(m, n, n2, maxrule, t, w, w1);
		xa = xc = a;
		xb = b;
		fx[0] = f.apply(b);
		fx[n] = f.apply(a);
		fev += 2;

		boolean gotoNEXT = true;
		while (true) {

			if (gotoNEXT) {
				n = 4;
				n2 = 2;
				m = mmax;
				order = 1;
				caution = xa < xc;
				h = xb - xa;
				k1 = h / (b - xa);
				if (k1 < 0.1) {
					k1 = 0.1;
				}
				h *= 0.5;
				j = 1;
				fx[n2] = f.apply(xa + h);
				++fev;
				fmin = fmax = fx[0];
				if (fmax < fx[n]) {
					fmax = fx[n];
				} else if (fmin > fx[n]) {
					fmin = fx[n];
				}
				if (fmax < fx[n2]) {
					fmax = fx[n2];
				} else if (fmin > fx[n2]) {
					fmin = fx[n2];
				}
			}

			// AGAIN
			for (i = 1; i <= n2 - 1; i += j) {
				fx[i] = f.apply(xa + (1.0 + t[i * m]) * h);
				++fev;
				if (fmax < fx[i]) {
					fmax = fx[i];
				} else if (fmin > fx[i]) {
					fmin = fx[i];
				}
				fx[n - i] = f.apply(xa + (1.0 - t[i * m]) * h);
				++fev;
				if (fmax < fx[n - i]) {
					fmax = fx[n - i];
				} else if (fmin > fx[n - i]) {
					fmin = fx[n - i];
				}
			}
			re = acc * Math.max(Math.abs(fmax), Math.abs(fmin));
			j = n == 4 ? 4 : 6;
			k = 0;
			c = Math.abs(cheb(-1.0, fx, n)) / n;
			for (i = 2; i <= j; i += 2) {
				x = i <= n2 ? -t[i * m] : t[(n - i) * m];
				cprev = c;
				c = Math.abs(cheb(x, fx, n)) / n2;
				if (c > re) {
					if (k < cprev / c) {
						k = cprev / c;
					}
				} else if (cprev > re) {
					k = 1;
				}
			}

			boolean gotoUPDATE = false;
			if (k > SIGMA[order - 1][4]) {
				if (n == 4) {
					e = 2.0 * h * (fmax - fmin);
					if (e < re || e <= k1 * eps) {
						iint = h * (fmax + fmin);
						gotoUPDATE = true;
					} else {
						if (h < hmin) {
							break;
						}
						xs[div] = xb;
						xb = xa + h;
						fs[div] = fx[0];
						fx[0] = fx[n2];
						fx[4] = fx[n];
						div++;
						gotoNEXT = true;
						continue;
					}
				}
			} else {
				if (n == 4) {
					cprev = c;
				} else if (cprev < re) {
					cprev = k * c;
				}
				e = h * cprev * ec[order - 1] * k * k * k;
				i = 0;
				while (k > SIGMA[order - 1][i]) {
					++i;
					e *= 2.0;
				}
				re *= h;
			}

			if (!gotoUPDATE) {
				iint = w1[n] * (fx[0] + fx[n]) + w[n - 3] * fx[n2];
				for (i = 1; i <= n2 - 1; ++i) {
					iint += w[n2 + i - 3] * (fx[i] + fx[n - i]);
				}
				iint *= h;
				if (n != 4) {
					c = Math.abs(iint - intprev);
					if (c > eprev) {
						caution = true;
						if (xc < xb) {
							xc = xb;
						}
					} else {
						caution = false;
					}
					if (k > SIGMA[order - 1][4] || caution) {
						e = c;
					}
					if (e > c) {
						e = c;
					}
				}

				// TEST
				if (e < re || e <= k1 * eps) {
				} else if (k > SIGMA[order - 1][0]) {
					e = 2.0 * h * (fmax - fmin);
					if (e < re || e <= k1 * eps) {
						iint = h * (fmax + fmin);
					} else {
						if (h < hmin) {
							break;
						}
						xs[div] = xb;
						xb = xa + h;
						fs[div] = fx[0];
						fx[0] = fx[n2];
						fx[4] = fx[n];
						++div;
						gotoNEXT = true;
						continue;
					}
				} else {
					for (i = n; i >= 1; --i) {
						fx[i << 1] = fx[i];
					}
					n2 = n;
					n <<= 1;
					m >>= 1;
					++order;
					eprev = e;
					intprev = iint;
					if (eprev < re) {
						eprev = re;
					}
					if (n > maxrule) {
						maxrule = quadrule(m, n, n2, maxrule, t, w, w1);
					}
					j = 2;
					gotoNEXT = false;
					continue;
				}
			}

			// UPDATE
			if (n != 4 || !(caution || (xa == a && div == 0))) {
				if (e < re) {
					e = re;
				}
				error += e;
				eps -= e;
				if (eps < 0.1 * error) {
					eps = 0.1 * error;
				}
				ans += iint;
				if (div == 0) {
					myFEvals += fev;
					return ans;
				}
				--div;
				xa = xb;
				xb = xs[div];
				fx[4] = fx[0];
				fx[0] = fs[div];
				gotoNEXT = true;
				continue;
			}

			// DOUBLE
			for (i = n; i >= 1; --i) {
				fx[i << 1] = fx[i];
			}
			n2 = n;
			n <<= 1;
			m >>= 1;
			++order;
			eprev = e;
			intprev = iint;
			if (eprev < re) {
				eprev = re;
			}
			if (n > maxrule) {
				maxrule = quadrule(m, n, n2, maxrule, t, w, w1);
			}
			j = 2;
			gotoNEXT = false;
		}
		myFEvals += fev;
		return Double.NaN;
	}

	private static int quadrule(final int m, final int n, final int n2, final int maxrule, final double[] t,
			final double[] w, final double[] w1) {
		for (int i = 1; i <= n2 - 1; i += 2) {
			t[i * m] = Math.cos(Constants.PI * i / n);
		}
		for (int i = (maxrule >> 1) + 2; i <= n; i += 2) {
			w1[i - 1] = 0.0;
			w1[i] = 1.0 / (i * i - 1.0);
		}
		for (int i = 1; i <= n2; ++i) {
			w[n2 + i - 3] = -4.0 * cheb(t[i * m], w1, n) / n;
		}
		return n;
	}

	private static double cheb(final double x, final double[] a, final int n) {
		double b0 = 0.5 * a[n];
		double b1 = 0.0;
		for (int r = n - 1; r >= 1; --r) {
			final double b2 = b1;
			b1 = b0;
			b0 = 2.0 * x * b1 - b2 + a[r];
		}
		return x * b0 - b1 + 0.5 * a[0];
	}
}
