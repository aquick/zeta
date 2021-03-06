package org.gmplib.test.zeta;

import java.util.ArrayList;

import org.gmplib.gmpjni.GMP;
import org.gmplib.gmpjni.MPFR;
import org.gmplib.gmpjni.GMP.GMPException;
import org.gmplib.gmpjni.GMP.mpz_t;
import org.gmplib.gmpjni.MPFR.MPFRException;
import org.gmplib.gmpjni.MPFR.mpfr_rnd_t;
import org.gmplib.gmpjni.MPFR.mpfr_t;

public class RiemannZetaBasic implements RiemannZeta {

    private StringBuffer zeroStr[];
    private int thetaEvals;
    private int zetaEvals;
    private int precision;
    private int digits; // precision in decimal digits
    private int numZeroes;
    
    private mpfr_t pi;
    private mpfr_t piby2;
    private mpfr_t piby8;
    private mpfr_t log2;
    private int numGramPoints = 0;
    private int gramBaseIndex = -1;
    private mpfr_t[] gramPoints; // gramPoints[i] is g(gramBaseIndex + i)
    private mpfr_t[] logs; // logs[i] is ln(i+1)
    private mpfr_t[] E; // coefficients for zeta_sum_part2 E[i] is e[i+1]/2^n
    private mpfr_t epsilon; // 1 / 2^(precision-2)
    private static final int[] Bnums = new int[] {1, 1, 1,  1,  1,  5,  691,  7, 3617, 43867, 174611, 854513}; // numerators of abs(B(2n))
    //                                      -  -  -   --  --  --  ----  -  ----  -----  ------  ------
    private static final int[] Bdens = new int[] {1, 6, 30, 42, 30, 66, 2730, 6, 510,  798,   330,    138   }; // denominators of abs(B(2n))
    private mpfr_t[] C; // coefficients for theta
    private mpfr_t[] tpowers; // used in theta

    // temporary variables used by mpfr_cmul, mpfr_cinvert
    private static mpfr_t c1;
    private static mpfr_t c2;
    // temporary variables used by theta, zeta_sum_part1, zeta_sum_part2, zeta_part3, computeNumberOfTerms, computeNumberOfZeroes
    private static mpfr_t t1;
    private static mpfr_t t2;
    private static mpfr_t t3;
    private static mpfr_t u;
    private static mpfr_t v;
    // temporary variables used by zeta
    private static mpfr_t x1;
    private static mpfr_t y1;
    private static mpfr_t sumx;
    private static mpfr_t sumy;
    // temporary variables used by Z
    private static mpfr_t t4;
    private static mpfr_t t5;
    private static mpfr_t thetax;
    private static mpfr_t thetay;
    private static mpfr_t zetax;
    private static mpfr_t zetay;

    /**
     * Constructor: initialize member variables, initialize constants.
     */
    public RiemannZetaBasic()
	throws MPFRException, Exception
    {	
	this.numZeroes = 0;
	thetaEvals = 0;
	zetaEvals = 0;
	epsilon = new mpfr_t();
	pi = new mpfr_t();
	piby2 = new mpfr_t();
	piby8 = new mpfr_t();
	log2 = new mpfr_t();
	MPFR.mpfr_const_pi(pi, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(piby2, pi, 2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(piby8, pi, 8, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_const_log2(log2, mpfr_rnd_t.MPFR_RNDN);
	c1 = new mpfr_t();
	c2 = new mpfr_t();
	t1 = new mpfr_t();
	t2 = new mpfr_t();
	t3 = new mpfr_t();
	t4 = new mpfr_t();
	t5 = new mpfr_t();
	x1 = new mpfr_t();
	y1 = new mpfr_t();
	sumx = new mpfr_t();
	sumy = new mpfr_t();
	u = new mpfr_t();
	v = new mpfr_t();
	thetax = new mpfr_t();
	thetay = new mpfr_t();
	zetax = new mpfr_t();
	zetay = new mpfr_t();
    }
    
    public void setPrecision(int precision)
	throws MPFRException, Exception
    {
        this.precision = precision;
        MPFR.mpfr_set_ui(epsilon, 1, mpfr_rnd_t.MPFR_RNDN);
        MPFR.mpfr_div_2si(epsilon, epsilon, precision - 2, mpfr_rnd_t.MPFR_RNDN);
        this.digits = (int)((double)(precision - 2)*0.301029995664 /* log10(2) */);
    }

    /**
     * Initialize array of logarithms.
     * logs[i] = ln(i+1) for i = 0,...,n-1
     */
    private void initLogs(int n)
        throws MPFRException
    {
	int i;
	logs = new mpfr_t[n];
	mpfr_t t = new mpfr_t();
	for (i = 0; i < n; i++) {
	    logs[i] = new mpfr_t();
	    MPFR.mpfr_set_ui(t, i + 1, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_log(logs[i], t, mpfr_rnd_t.MPFR_RNDN);
	}
    }

    /**
     * Initialize basic functions.  In this case logarithms.
     */
    public void initBasicFunctions(int n, mpfr_t tmin, mpfr_t tmax)
	throws MPFRException
    {
	initLogs(2*n);
    }

    /**
     * Pre-compute coefficients used by zeta_sum_part2.
     * E[i] is e[i+1]/2^n, where e[k]=sum (n choose j) j=k,...,n
     */
    public void initCoefficients(int n)
        throws MPFRException, GMPException
    {
	int k;
	int j = n - 1;
	mpz_t e = new mpz_t();
	mpz_t f = new mpz_t();
	E = new mpfr_t[n];
	GMP.mpz_set_si(f, 1);
	GMP.mpz_set(e, f);
	for (k = 2*n; k >= n + 1; k--) {
	    E[j] = new mpfr_t();
	    MPFR.mpfr_set_z(u, e, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_div_2si(u, u, n, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_set(E[j], u, mpfr_rnd_t.MPFR_RNDN);
	    j--;
            GMP.mpz_mul_si(f, f, k - n);
	    GMP.mpz_divexact_ui(f, f, (long)(2*n - k + 1));
	    GMP.mpz_add(e, e, f);
	}
    }

    /**
     * Pre-compute coefficients used by theta.
     * C[n]=(1/2 - 1/2^2*n)*(1/((2*n-1)*2*n)*abs(B[2*n]) where B[i] is the i-th Bernoulli number.
     */
    public void initThetaCoefficients()
        throws MPFRException
    {
	mpfr_t c = new mpfr_t();
	mpfr_t d = new mpfr_t();
	int i;

	C = new mpfr_t[Bnums.length - 1];
	tpowers = new mpfr_t[Bnums.length - 1];

	MPFR.mpfr_set_si(c, 1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_2si(c, c, 1, mpfr_rnd_t.MPFR_RNDN);
	for (i = 1; i < Bnums.length; i++) {
	    C[i-1] = new mpfr_t();
	    MPFR.mpfr_set_si(d, 1, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_div_2si(d, d, 2*i, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_sub(d, c, d, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul_si(d, d, Bnums[i], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_div_si(d, d, Bdens[i], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_div_si(d, d, 2*i - 1, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_div_si(d, d, 2*i, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_set(C[i-1], d, mpfr_rnd_t.MPFR_RNDN);
	    tpowers[i-1] = new mpfr_t();
	}
    }

    /**
     * Multiply complex numbers x1 + i*y1 and x2 + i*y2.
     */
    private static void mpfr_cmul(mpfr_t rx, mpfr_t ry, mpfr_t x1, mpfr_t y1, mpfr_t x2, mpfr_t y2, mpfr_rnd_t rnd)
        throws MPFRException
    {
	MPFR.mpfr_mul(rx, x1, x2, rnd);
	MPFR.mpfr_mul(c1, y1, y2, rnd);
	MPFR.mpfr_sub(rx, rx, c1, rnd);
	MPFR.mpfr_mul(ry, x1, y2, rnd);
	MPFR.mpfr_mul(c1, x2, y1, rnd);
	MPFR.mpfr_add(ry, ry, c1, rnd);
    }

    /**
     * Invert a complex number x + i*y.
     */
    private static void mpfr_cinvert(mpfr_t rx, mpfr_t ry, mpfr_t x, mpfr_t y, mpfr_rnd_t rnd)
        throws MPFRException
    {
	MPFR.mpfr_sqr(c1, x, rnd);
	MPFR.mpfr_sqr(c2, y, rnd);
	MPFR.mpfr_add(c1, c1, c2, rnd);
	MPFR.mpfr_div(rx, x, c1, rnd);
	MPFR.mpfr_div(ry, y, c1, rnd);
	MPFR.mpfr_neg(ry, ry, rnd);
    }
    
    /**
     * Compute high-order terms of the Riemann-Siegel theta function for t >> 1.
     * Used by function theta.  Does not initialize r.
     */
    private void theta_high_order_terms(mpfr_t r, mpfr_t t)
	        throws MPFRException
    {
	MPFR.mpfr_sub(r, r, piby8, mpfr_rnd_t.MPFR_RNDN);
	
	MPFR.mpfr_set(v,  t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_2si(v, v, 1, mpfr_rnd_t.MPFR_RNDN);
        MPFR.mpfr_sub(r, r, v, mpfr_rnd_t.MPFR_RNDN);
        
	MPFR.mpfr_set(u, t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_2si(u, u, 1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(u, u, pi, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(u, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(v, v, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(r, r, v, mpfr_rnd_t.MPFR_RNDN);	
    }

    /**
     * Approximate the Riemann-Siegel theta function for t >> 1.
     * Reference: Dissertation: NEUE HERLEITUNG UND EXPLIZITE RESTABSCHATZUNG DER RIEMANN-SIEGEL-FORMEL,
     * Wolfgang Gabcke, Gottingen 1979, Page 2.
     */
    private void theta(mpfr_t r, mpfr_t t, int nterms)
        throws MPFRException
    {
	int i;
	
	thetaEvals++;
	MPFR.mpfr_sqr(t2, t, mpfr_rnd_t.MPFR_RNDN);	
	for (i = 0; i < nterms; i++) {
	    if (i == 0) {
		MPFR.mpfr_set_si(tpowers[i], 1, mpfr_rnd_t.MPFR_RNDN);
		MPFR.mpfr_div(tpowers[i], tpowers[i], t, mpfr_rnd_t.MPFR_RNDN);
	    } else {
		MPFR.mpfr_div(tpowers[i], tpowers[i-1], t2, mpfr_rnd_t.MPFR_RNDN);
	    }
	}
	MPFR.mpfr_set_ui(r, 0, mpfr_rnd_t.MPFR_RNDN);
	for (i = nterms - 1; i >= 0; i--) {
	    MPFR.mpfr_mul(v, C[i], tpowers[i], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_add(r, r, v, mpfr_rnd_t.MPFR_RNDN);
	}

	theta_high_order_terms(r, t);
    }

    /**
     * Initialize a list of Gram points between lb and ub.
     * Gram points are found by binary search for solutions to theta(t) = n*pi.
     */
    public void initGramPoints(int lb, int ub, int nt)
        throws MPFRException
    {
	ArrayList<mpfr_t> points = new ArrayList<mpfr_t>();
	int i = 0;
	mpfr_t s = new mpfr_t();
	mpfr_t a = new mpfr_t();
	mpfr_t b = new mpfr_t();
	mpfr_t d = new mpfr_t();
	mpfr_t mid = new mpfr_t();
	mpfr_t r = new mpfr_t();
	mpfr_t p;
        MPFR.mpfr_neg(s, pi, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_si(a, lb, mpfr_rnd_t.MPFR_RNDN);
	theta(r, a, nt);
	for (;;) {
    	    if (MPFR.mpfr_cmp(r, s) < 0) {
    	        break;
    	    }
            MPFR.mpfr_add(s, s, pi, mpfr_rnd_t.MPFR_RNDN);
            i++;
    	}
        gramBaseIndex = i - 1;	    
	for (;;) {
	    p = new mpfr_t();
            for (;;) {
        	MPFR.mpfr_add_si(b, a, 1, mpfr_rnd_t.MPFR_RNDN);
        	theta(r, b, nt);
        	if (MPFR.mpfr_cmp(r, s) > 0) {
        	    break;
        	}
        	MPFR.mpfr_set(a, b, mpfr_rnd_t.MPFR_RNDN);
            }
            // now we have theta(a) < s < theta(b)
            for (;;) {
        	MPFR.mpfr_add(mid, a, b, mpfr_rnd_t.MPFR_RNDN);
        	MPFR.mpfr_div_si(mid, mid, 2, mpfr_rnd_t.MPFR_RNDN);
        	theta(r, mid, nt);
        	if (MPFR.mpfr_cmp(r, s) < 0) {
        	    MPFR.mpfr_set(a, mid, mpfr_rnd_t.MPFR_RNDD);
        	} else {
        	    MPFR.mpfr_set(b, mid, mpfr_rnd_t.MPFR_RNDU);
        	}
        	MPFR.mpfr_div(d, a, b, mpfr_rnd_t.MPFR_RNDN);
        	MPFR.mpfr_d_sub(d, 1.0, d, mpfr_rnd_t.MPFR_RNDN);
        	if (MPFR.mpfr_cmp_d(d, 0.00001) < 0) { // 1E-5
        	    break;
        	}
            }
            MPFR.mpfr_set(p, a, mpfr_rnd_t.MPFR_RNDN);
            if (MPFR.mpfr_cmp_si(p, ub) > 0) {
        	break;
            }
            points.add(p);
            MPFR.mpfr_add(s, s, pi, mpfr_rnd_t.MPFR_RNDN);
            i++;
	}
	gramPoints = points.toArray(new mpfr_t[0]);
	numGramPoints = points.size();
        zeroStr = new StringBuffer[numGramPoints + 1];
        for (i = 0; i <= numGramPoints; i++) {
            zeroStr[i] = new StringBuffer();
        }
    }
    
    /**
     * Returns the number of Gram points computed by initGramPoints.
     */
    public int getNumberOfGramPoints()
    {
	return this.numGramPoints;
    }

    /**
     * Returns the index of the first Gram point computed by initGramPoints.
     */
    public int getGramBaseIndex()
    {
	return this.gramBaseIndex;
    }

    /**
     * Returns the i-th Gram point computed by initGramPoints.
     */
    public mpfr_t getGramPoint(int i)
	throws MPFRException
    {
	return this.gramPoints[i];
    }

    /**
     * Compute alternating sum (-1)^(k-1)/k^s, k=1,...,n
     * Used by function zeta.
     */
    private void zeta_sum_part1(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy, int n)
        throws MPFRException
    {
	int k;
	MPFR.mpfr_set_si(rx,  0, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_si(ry,  0, mpfr_rnd_t.MPFR_RNDN);
	for (k = 1; k <= n; k++) {
	    MPFR.mpfr_mul(u, sy, logs[k-1], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_neg(u, u, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(v, sx, logs[k-1], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_neg(v, v, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_sin_cos(t1, t2, u, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_exp(t3, v, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(t1, t1, t3, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(t2, t2, t3, mpfr_rnd_t.MPFR_RNDN);
	    if (k % 2 == 0) {
		MPFR.mpfr_sub(rx, rx, t2, mpfr_rnd_t.MPFR_RNDN);
		MPFR.mpfr_sub(ry, ry, t1, mpfr_rnd_t.MPFR_RNDN);
	    } else {
		MPFR.mpfr_add(rx, rx, t2, mpfr_rnd_t.MPFR_RNDN);
		MPFR.mpfr_add(ry, ry, t1, mpfr_rnd_t.MPFR_RNDN);
	    }
	}
    }
    
    /**
     * Compute alternating sum (-1)^(k-1)*E[k-n]/k^s, k=n+1,...,2*n
     * Used by function zeta.
     */
    private void zeta_sum_part2(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy, int n)
        throws MPFRException
    {
	int k;
	int j = n - 1;
	MPFR.mpfr_set_si(rx,  0, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_si(ry,  0, mpfr_rnd_t.MPFR_RNDN);
	for (k = 2*n; k >= n + 1; k--) {
	    MPFR.mpfr_mul(u, sy, logs[k-1], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_neg(u, u, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(v, sx, logs[k-1], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_neg(v, v, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_sin_cos(t1, t2, u, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_exp(t3, v, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(t1, t1, t3, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(t2, t2, t3, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(t1, t1, E[j], mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_mul(t2, t2, E[j], mpfr_rnd_t.MPFR_RNDN);
	    j--;
	    if (k % 2 == 0) {
		MPFR.mpfr_sub(rx, rx, t2, mpfr_rnd_t.MPFR_RNDN);
		MPFR.mpfr_sub(ry, ry, t1, mpfr_rnd_t.MPFR_RNDN);
	    } else {
		MPFR.mpfr_add(rx, rx, t2, mpfr_rnd_t.MPFR_RNDN);
		MPFR.mpfr_add(ry, ry, t1, mpfr_rnd_t.MPFR_RNDN);
	    }
	}
    }

    /**
     * Compute 1/(1 - 2^(1 - s)) where s = sx + i*sy
     * Used by function zeta.
     */
    private void zeta_part3(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy)
        throws MPFRException
    {
        MPFR.mpfr_mul(u, sy, log2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_neg(u, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_ui_sub(v, 1, sx, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(v, v, log2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_sin_cos(t1, t2, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_exp(t3, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(t1, t1, t3, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(t2, t2, t3, mpfr_rnd_t.MPFR_RNDN);
	// t2 = Re(2^(1-s))
	// t1 = Im(2^(1-s))
	MPFR.mpfr_neg(t1, t1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_neg(t2, t2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add_ui(t2, t2, 1, mpfr_rnd_t.MPFR_RNDN);
	mpfr_cinvert(rx, ry, t2, t1, mpfr_rnd_t.MPFR_RNDN);
    }
    
    /**
     * Compute zeta(s) where s = sx + i*sy and sx > 0.
     * Reference: Numerical evaluation of the Riemann Zeta-function, Xavier Gourdon and Pascal Sebah, July 23, 2003.
     * Alternating series method of Proposition 2.
     * @see http://numbers.computation.free.fr/Constants/constants.html
     */
    private void zeta(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy, int n)
        throws MPFRException
    {
	zeta_sum_part1(sumx, sumy, sx, sy, n);
	zeta_sum_part2(x1, y1, sx, sy, n);
	MPFR.mpfr_add(sumx, sumx, x1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(sumy, sumy, y1, mpfr_rnd_t.MPFR_RNDN);
	zeta_part3(x1, y1, sx, sy);
	mpfr_cmul(rx, ry, sumx, sumy, x1, y1, mpfr_rnd_t.MPFR_RNDN);
    }
    
    /**
     * Compute Z(t) = exp(i*theta(t)) * zeta(0.5 + i*t) for t > 0
     */
    public void Z(mpfr_t r, mpfr_t t, int n, int nt)
        throws MPFRException, Exception
    {
	String str;
        GMP.MutableInteger exp = new GMP.MutableInteger(0);

	zetaEvals++;
        theta(t4, t, nt);
	MPFR.mpfr_sin_cos(thetay, thetax, t4, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_d(t4, 0.5, mpfr_rnd_t.MPFR_RNDN);
	zeta(zetax, zetay, t4, t, n);
	mpfr_cmul(t4, t5, thetax, thetay, zetax, zetay, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_abs(t5, t5, mpfr_rnd_t.MPFR_RNDN);
	if (MPFR.mpfr_cmp_d(t5, 0.01) > 0) {
	    str = MPFR.mpfr_get_str(exp, t5, 10, 8, mpfr_rnd_t.MPFR_RNDN);
	    throw new Exception("Z: abs(Im(Z))=0." + str + "E" + Integer.toString(exp.value));
	}
	MPFR.mpfr_set(r, t4, mpfr_rnd_t.MPFR_RNDN);
    }
    
    /**
     * Append a base 10 string representation of an mpfr_t to StringBuffer result.  The number of
     * digits is based on the desired precision.
     */
    private void appendMPFRvalue(StringBuffer sb, mpfr_t value)
        throws MPFRException
    {
        GMP.MutableInteger exp = new GMP.MutableInteger(0);
	String str = MPFR.mpfr_get_str(exp, value, 10, this.digits, mpfr_rnd_t.MPFR_RNDN);
	if (str.startsWith("-")) {
	    sb.append("-0.");
	    sb.append(str.substring(1));
	} else {
	    sb.append("0.");
	    sb.append(str);
	}
	sb.append("E");
	sb.append(Integer.toString(exp.value));	
    }
    
    /**
     * Append a base 10 string representation of an mpfr_t to StringBuffer result.  The number of
     * digits is given by ndigits.
     */
    private void appendMPFRvalue(StringBuffer sb, mpfr_t value, int ndigits)
        throws MPFRException
    {
        GMP.MutableInteger exp = new GMP.MutableInteger(0);
	String str = MPFR.mpfr_get_str(exp, value, 10, ndigits, mpfr_rnd_t.MPFR_RNDN);
	if (str.startsWith("-")) {
	    sb.append("-0.");
	    sb.append(str.substring(1));
	} else {
	    sb.append("0.");
	    sb.append(str);
	}
	sb.append("E");
	sb.append(Integer.toString(exp.value));	
    }

    /**
     * Compute a zero of zeta(0.5 + i*t) with a0 < t < b0.
     * If there is a sign change of Z in the interval then linear interpolation is used to locate a zero.
     * @return true iff a zero was found in the interval.
     */
    public boolean computeZetaZero(int i, mpfr_t a0, mpfr_t b0, int n, int nt)
        throws MPFRException, Exception
    {
	mpfr_t a = new mpfr_t();
	mpfr_t b = new mpfr_t();
	mpfr_t d = new mpfr_t();
	mpfr_t mid = new mpfr_t();
	mpfr_t zmid = new mpfr_t();
	mpfr_t za = new mpfr_t();
	mpfr_t zb = new mpfr_t();
	boolean up = false; // true if Z(a) < 0 < Z(b); false if Z(a) > 0 > Z(b) 

	MPFR.mpfr_set(a, a0, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set(b, b0, mpfr_rnd_t.MPFR_RNDN);
	Z(za, a, n, nt);
	Z(zb, b, n, nt);
	if (MPFR.mpfr_cmp_ui(za, 0) < 0 && MPFR.mpfr_cmp_ui(zb, 0) > 0) {
	    up = true;
	} else if (MPFR.mpfr_cmp_ui(za, 0) > 0 && MPFR.mpfr_cmp_ui(zb, 0) < 0) {
	    up = false;
	} else {
	    this.zeroStr[i].append("[");
	    appendMPFRvalue(this.zeroStr[i], a);
	    this.zeroStr[i].append(", ");
	    appendMPFRvalue(this.zeroStr[i], b);
	    this.zeroStr[i].append("]");
	    return false;
	}
	for (;;) {
	    // set mid to the zero of the line between (a, Z(a)) and (b, Z(b))
            MPFR.mpfr_sub(mid, b, a, mpfr_rnd_t.MPFR_RNDN);
            MPFR.mpfr_div(d, zb, za, mpfr_rnd_t.MPFR_RNDN);
            MPFR.mpfr_d_sub(d, 1.0, d, mpfr_rnd_t.MPFR_RNDN);
            if (MPFR.mpfr_cmp_d(d, 1.0) <= 0) {
        	throw new Exception("computeZetaZero: 1 - zb/za <= 1");        	
            }
            MPFR.mpfr_div(mid, mid, d, mpfr_rnd_t.MPFR_RNDN);
            MPFR.mpfr_add(mid, mid, a, mpfr_rnd_t.MPFR_RNDN);
            // if a == mid or b == mid, we take mid as the average of a and b
            if (MPFR.mpfr_cmp(a, mid) == 0 || MPFR.mpfr_cmp(b, mid) == 0) {
        	MPFR.mpfr_add(mid, a, b, mpfr_rnd_t.MPFR_RNDN);
        	MPFR.mpfr_div_2si(mid, mid, 1, mpfr_rnd_t.MPFR_RNDN);
            }
            Z(zmid, mid, n, nt);
            if (up == (MPFR.mpfr_cmp_ui(zmid, 0) < 0)) {
        	MPFR.mpfr_set(a, mid, mpfr_rnd_t.MPFR_RNDD);
        	MPFR.mpfr_set(za, zmid, (up ? mpfr_rnd_t.MPFR_RNDD : mpfr_rnd_t.MPFR_RNDU));
            } else {
        	MPFR.mpfr_set(b, mid, mpfr_rnd_t.MPFR_RNDU);
        	MPFR.mpfr_set(zb, zmid, (up ? mpfr_rnd_t.MPFR_RNDU : mpfr_rnd_t.MPFR_RNDD));
            }
            MPFR.mpfr_div(d, a, b, mpfr_rnd_t.MPFR_RNDN);
            MPFR.mpfr_d_sub(d, 1.0, d, mpfr_rnd_t.MPFR_RNDN);
            // break if the relative difference between a and b is within epsilon
            if (MPFR.mpfr_cmp(d, epsilon) < 0) {
        	break;
            }
        }
	appendMPFRvalue(this.zeroStr[i], a);
	this.numZeroes++;
	return true;
    }

    /**
     * Returns the total number of zeroes computed by computeZetaZero.
     */
    public int getNumberOfZeroes()
    {
	return this.numZeroes;
    }

    /**
     * Compute the number of terms needed in zeta sums for the desired precision.
     */
    public int computeNumberOfTerms(mpfr_t tmin, mpfr_t tmax)
        throws MPFRException
    {
	MPFR.mpfr_mul_ui(t1, tmax, 2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add_ui(t1, t1, 1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(u, t1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(t1, tmax, piby2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(u, u, t1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul_ui(t1, log2, this.precision - 2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(u, u, t1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_sub_d(u, u, 0.23054022971699147, mpfr_rnd_t.MPFR_RNDN); // ln(sqrt(3 - 2*sqrt(2))) = 0.23054022971699147077486983704068
	MPFR.mpfr_mul_ui(t1, log2, 3, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(u, u, t1, mpfr_rnd_t.MPFR_RNDN);
	return MPFR.mpfr_get_si(u, mpfr_rnd_t.MPFR_RNDD);
    }

    /**
     * Compute the number of terms needed in the theta function for the desired precision.
     */
    public int computeNumberOfThetaTerms(mpfr_t tmin)
        throws MPFRException
    {
	int i;
	int n = 0;
	mpfr_t tp = new mpfr_t();
	mpfr_t x = new mpfr_t();
	mpfr_t y = new mpfr_t();
	mpfr_t d = new mpfr_t();

	MPFR.mpfr_set_si(y, 0, mpfr_rnd_t.MPFR_RNDN);
	theta_high_order_terms(y, tmin);
	MPFR.mpfr_sqr(t2, tmin, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_si(tp, 1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(tp, tp, tmin, mpfr_rnd_t.MPFR_RNDN);
	for (i = 0; i < C.length; i++) {
	    MPFR.mpfr_mul(x, C[i], tp, mpfr_rnd_t.MPFR_RNDN);
    	    MPFR.mpfr_div(d, x, y, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_abs(d, d, mpfr_rnd_t.MPFR_RNDN);
	    if (MPFR.mpfr_cmp(d, epsilon) < 0) break;
	    MPFR.mpfr_div(tp, tp, t2, mpfr_rnd_t.MPFR_RNDN);
	    n++;
	}
	return n;
    }

    /**
     * Compute an estimate of N(T), the number of zeroes of zeta(0.5 + i*y) for i in [0,t]
     * from Backlund.
     * Reference: Riemann's Zeta Function, H. M. Edwards, Dover Publications Inc. 2001.
     */
    public long computeNumberOfZeroes(mpfr_t error, mpfr_t t)
        throws MPFRException, Exception
    {
	long n;
	mpfr_t nt = new mpfr_t();
	
	// hard-code the first few estimates
	if (MPFR.mpfr_cmp_d(t, 14.134725142) < 0) {
	    n = 0;
	    MPFR.mpfr_set_si(error, 0, mpfr_rnd_t.MPFR_RNDN);
	    return n;
	} else if (MPFR.mpfr_cmp_d(t, 21.022039639) < 0) {
	    n = 1;
	    MPFR.mpfr_set_si(error, 0, mpfr_rnd_t.MPFR_RNDN);
	    return n;
	} else if (MPFR.mpfr_cmp_d(t, 25.010857580) < 0) {
	    n = 2;
	    MPFR.mpfr_set_si(error, 0, mpfr_rnd_t.MPFR_RNDN);
	    return n;
	}
	// find the largest Gram point less than the bound and, if the index is 126 or less,
	// use the well known property of Gram's Law
	int i = 0;
	while (i < numGramPoints) {
	    if (MPFR.mpfr_cmp(gramPoints[i], t) >= 0) {
		break;
	    }
	    i++;
	}
	if (gramBaseIndex + i <= 126) {
	    n = gramBaseIndex + i;
	    MPFR.mpfr_set_si(error, 1, mpfr_rnd_t.MPFR_RNDN);
	    return n;
	}
	// otherwise use Backlund's estimate
	MPFR.mpfr_set(v,  t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(v, v, 2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(v, v, pi, mpfr_rnd_t.MPFR_RNDN);
        
	MPFR.mpfr_set(u, t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(u, u, 2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(u, u, pi, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(u, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(u, u, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set(nt, u, mpfr_rnd_t.MPFR_RNDN);
	
	MPFR.mpfr_sub(nt, nt, v, mpfr_rnd_t.MPFR_RNDN);
	
	MPFR.mpfr_add_d(nt, nt, 0.875, mpfr_rnd_t.MPFR_RNDN);

	n = MPFR.mpfr_get_ui(nt, mpfr_rnd_t.MPFR_RNDN);
	if (MPFR.mpfr_erangeflag_p() != 0) {
	    throw new Exception("Overflow in computeNumberOfZeroes");
	}

	MPFR.mpfr_set(v,  t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(v, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul_d(error, v, 0.137, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(v, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul_d(v, v, 0.443, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(error, error, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add_d(error, error, 4.35, mpfr_rnd_t.MPFR_RNDN);
	return n;
    }

    /**
     * Returns a string result for the i-th call to computeZetaZero.
     */
    public String getZeroStr(int i)
    {
	return this.zeroStr[i].toString();
    }

    /**
     * Returns the total number of evaluations of Z.
     */
    public int getNumberOfZEvals()
    {
	return this.zetaEvals;
    }
    
    /**
     * Returns the total number of evaluations of theta.
     */
    public int getNumberOfThetaEvals()
    {
	return this.thetaEvals;
    }

}
