package org.gmplib.test.zeta;

import android.os.AsyncTask;
import android.util.Log;

import java.util.ArrayList;

import org.gmplib.gmpjni.GMP;
import org.gmplib.gmpjni.GMP.mpz_t;
import org.gmplib.gmpjni.GMP.GMPException;
import org.gmplib.gmpjni.MPFR;
import org.gmplib.gmpjni.MPFR.mpfr_t;
import org.gmplib.gmpjni.MPFR.mpfr_rnd_t;
import org.gmplib.gmpjni.MPFR.MPFRException;

public class Zeta_Task extends AsyncTask<Integer, Integer, Integer>
{

    private static final String TAG = "Zeta_Task";
    private UI uinterface;
    private StringBuffer result;
    private StringBuffer zeroStr[];
    private static int thetaEvals;
    private static int zetaEvals;
    private int precision;
    private int digits; // precision in decimal digits
    private int numZeroes;
    private int lowerBound;
    private int upperBound;
    
    private static mpfr_t pi;
    private static mpfr_t piby2;
    private static mpfr_t piby8;
    private static mpfr_t log2;
    private static int numGramPoints = 0;
    private static int gramBaseIndex = -1;
    private static mpfr_t[] gramPoints; // gramPoints[i] is g(gramBaseIndex + i)
    private static mpfr_t[] logs; // logs[i] is ln(i+1)
    private static mpfr_t[] E; // coefficients for zeta_sum_part2 E[i] is e[i+1]/2^n
    private static mpfr_t epsilon; // 1 / 2^(precision-2)

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

    public Zeta_Task(UI ui)
        throws MPFRException, Exception
    {
	this.uinterface = ui;
	if (ui.getPrecision() == 0) {
	    this.precision = MPFR.mpfr_get_default_prec();
	} else {
	    this.precision = ui.getPrecision();
	}
	if (this.precision < 4) {
	    throw new Exception("precision too low (must be at least 4)");
	}
	MPFR.mpfr_set_default_prec(this.precision);
	this.result = new StringBuffer();
	this.numZeroes = 0;
	this.lowerBound = 0;
	this.upperBound = 0;
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

    private static void initLogs(int n)
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
    
    private static void initCoefficients(int n)
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

    // approximate Riemann-Siegel theta function for t >> 1
    private static void theta(mpfr_t r, mpfr_t t)
        throws MPFRException
    {
	thetaEvals++;
	MPFR.mpfr_sqr(t2, t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_sqr(u, t2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(u, u, t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_ui(v, 31, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(v, v, 80640, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(v, v, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set(r, v, mpfr_rnd_t.MPFR_RNDN);
	
	MPFR.mpfr_set(u, t2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(u, u, t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_ui(v, 7, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(v, v, 5760, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(v, v, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(r, r, v, mpfr_rnd_t.MPFR_RNDN);

	MPFR.mpfr_set(u, t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_ui(v, 1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(v, v, 48, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(v, v, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(r, r, v, mpfr_rnd_t.MPFR_RNDN);
	
	MPFR.mpfr_sub(r, r, piby8, mpfr_rnd_t.MPFR_RNDN);
	
	MPFR.mpfr_set(v,  t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(v, v, 2, mpfr_rnd_t.MPFR_RNDN);
        MPFR.mpfr_sub(r, r, v, mpfr_rnd_t.MPFR_RNDN);
        
	MPFR.mpfr_set(u, t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div_ui(u, u, 2, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_div(u, u, pi, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(u, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul(v, v, u, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(r, r, v, mpfr_rnd_t.MPFR_RNDN);
    }
    
    private static void initGramPoints(int lb, int ub)
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
	theta(r, a);
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
        	theta(r, b);
        	if (MPFR.mpfr_cmp(r, s) > 0) {
        	    break;
        	}
        	MPFR.mpfr_set(a, b, mpfr_rnd_t.MPFR_RNDN);
            }
            // now we have theta(a) < s < theta(b)
            for (;;) {
        	MPFR.mpfr_add(mid, a, b, mpfr_rnd_t.MPFR_RNDN);
        	MPFR.mpfr_div_si(mid, mid, 2, mpfr_rnd_t.MPFR_RNDN);
        	theta(r, mid);
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
    }
    
    private static void zeta_sum_part1(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy, int n)
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
    
    private static void zeta_sum_part2(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy, int n)
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

    // compute 1/(1 - 2^(1 - s)) where s = sx + i*sy
    private static void zeta_part3(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy)
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
    
    // compute zeta(s) where s = sx + i*sy and sx > 0
    private static void zeta(mpfr_t rx, mpfr_t ry, mpfr_t sx, mpfr_t sy, int n)
        throws MPFRException
    {
	zetaEvals++;
	zeta_sum_part1(sumx, sumy, sx, sy, n);
	zeta_sum_part2(x1, y1, sx, sy, n);
	MPFR.mpfr_add(sumx, sumx, x1, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(sumy, sumy, y1, mpfr_rnd_t.MPFR_RNDN);
	zeta_part3(x1, y1, sx, sy);
	mpfr_cmul(rx, ry, sumx, sumy, x1, y1, mpfr_rnd_t.MPFR_RNDN);
    }
    
    // compute Z(t) = exp(i*theta(t)) * zeta(0.5 + i*t) for t > 0
    private static void Z(mpfr_t r, mpfr_t t, int n)
        throws MPFRException
    {
	String str;
        GMP.MutableInteger exp = new GMP.MutableInteger(0);

        theta(t4, t);
	MPFR.mpfr_sin_cos(thetay, thetax, t4, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_set_d(t4, 0.5, mpfr_rnd_t.MPFR_RNDN);
	zeta(zetax, zetay, t4, t, n);
	mpfr_cmul(t4, t5, thetax, thetay, zetax, zetay, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_abs(t5, t5, mpfr_rnd_t.MPFR_RNDN);
	if (MPFR.mpfr_cmp_d(t5, 0.01) > 0) {
	    str = MPFR.mpfr_get_str(exp, t5, 10, 8, mpfr_rnd_t.MPFR_RNDN);
	    throw new MPFRException(MPFRException.INTERNAL_ERROR, "abs(Im(Z))=0." + str + "E" + Integer.toString(exp.value));
	}
	MPFR.mpfr_set(r, t4, mpfr_rnd_t.MPFR_RNDN);
    }
    
    private void dumpZeta(double x, double y, int n)
        throws MPFRException
    {
        mpfr_t zx = new mpfr_t();
        mpfr_t zy = new mpfr_t();
        mpfr_t sx = new mpfr_t();
        mpfr_t sy = new mpfr_t();

        MPFR.mpfr_set_d(sx, 0.5, mpfr_rnd_t.MPFR_RNDN);
        MPFR.mpfr_set_d(sy, y, mpfr_rnd_t.MPFR_RNDN);
        zeta(zx, zy, sx, sy, n);
        appendMPFRvalue(this.result, zx);
	result.append(" + (");
        appendMPFRvalue(this.result, zy);
	result.append(")i");
	
	result.append("\t[Z]");
	Z(zy, sy, n);
        appendMPFRvalue(this.result, zy);
	result.append("\n");
    }

    // compute a zero of zeta(0.5 + i*t) with a0 < t < b0
    private boolean computeZetaZero(int i, mpfr_t a0, mpfr_t b0, int n)
        throws MPFRException
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
	Z(za, a, n);
	Z(zb, b, n);
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
        	throw new MPFRException(MPFRException.INTERNAL_ERROR, "1 - zb/za <= 1");        	
            }
            MPFR.mpfr_div(mid, mid, d, mpfr_rnd_t.MPFR_RNDN);
            MPFR.mpfr_add(mid, mid, a, mpfr_rnd_t.MPFR_RNDN);
            // if a == mid or b == mid, we take mid as the average of a and b
            if (MPFR.mpfr_cmp(a, mid) == 0 || MPFR.mpfr_cmp(b, mid) == 0) {
        	MPFR.mpfr_add(mid, a, b, mpfr_rnd_t.MPFR_RNDN);
        	MPFR.mpfr_div_2si(mid, mid, 1, mpfr_rnd_t.MPFR_RNDN);
            }
            Z(zmid, mid, n);
            if (up == (MPFR.mpfr_cmp_ui(zmid, 0) < 0)) {
        	MPFR.mpfr_set(a, mid, mpfr_rnd_t.MPFR_RNDD);
        	MPFR.mpfr_set(za, zmid, mpfr_rnd_t.MPFR_RNDN);
            } else {
        	MPFR.mpfr_set(b, mid, mpfr_rnd_t.MPFR_RNDU);
        	MPFR.mpfr_set(zb, zmid, mpfr_rnd_t.MPFR_RNDN);
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
    
    // Compute number of terms needed in zeta sums
    private int computeNumberOfTerms(mpfr_t tmax)
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

    // Compute an estimate of N(T) from Backlund
    private static void computeNumberOfZeroes(mpfr_t nt, mpfr_t error, mpfr_t t, boolean lb)
        throws MPFRException
    {
	if (lb && gramBaseIndex <= 126) {
	    MPFR.mpfr_set_si(nt, gramBaseIndex + 1, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_set_si(error, 1, mpfr_rnd_t.MPFR_RNDN);
	    return;
	}
	if (!lb && gramBaseIndex + numGramPoints - 1 <= 126) {
	    MPFR.mpfr_set_si(nt, gramBaseIndex + numGramPoints, mpfr_rnd_t.MPFR_RNDN);
	    MPFR.mpfr_set_si(error, 1, mpfr_rnd_t.MPFR_RNDN);
	    return;
	}
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

	MPFR.mpfr_set(v,  t, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(v, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul_d(error, v, 0.137, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_log(v, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_mul_d(v, v, 0.443, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add(error, error, v, mpfr_rnd_t.MPFR_RNDN);
	MPFR.mpfr_add_d(error, error, 4.35, mpfr_rnd_t.MPFR_RNDN);
    }

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
	    
    protected Integer doInBackground(Integer... params)
    {
	int i;
        int rc = -1;
        int nterms = 20;
        
        if (params.length > 0) {
            this.lowerBound = params[0].intValue();
        }
        if (params.length > 1) {
            this.upperBound = params[1].intValue();
        }
        try {
            if (this.lowerBound < 10) {
                throw new Exception("lower bound must be at least 10");
            }
            if (this.lowerBound >= this.upperBound) {
                throw new Exception("lower bound must be less than upper bound");
            }
            mpfr_t x = new mpfr_t();
            mpfr_t nzeroes = new mpfr_t();
            mpfr_t error = new mpfr_t();
            MPFR.mpfr_set_ui(epsilon, 1, mpfr_rnd_t.MPFR_RNDN);
            MPFR.mpfr_div_2si(epsilon, epsilon, precision - 2, mpfr_rnd_t.MPFR_RNDN);
            this.digits = (int)((double)(precision - 2)*0.301029995664 /* log10(2) */);
            initGramPoints(this.lowerBound, this.upperBound);
            result.append("Gram Points\n");
            for (i = 0; i < numGramPoints; i++) {
        	result.append("    g(" + Integer.toString(gramBaseIndex + i) + ")=");
        	appendMPFRvalue(this.result, gramPoints[i], 5);
        	result.append("\n");
            }
            // set nterms for the zero x + iy with maximum y
            // (if nterms is changed we must call initCoefficients again)
            MPFR.mpfr_set_si(x, this.upperBound, mpfr_rnd_t.MPFR_RNDN);
            nterms = computeNumberOfTerms(x);
            initLogs(2*nterms);
            initCoefficients(nterms);
            result.append("using 2*" + nterms + " terms in zeta sums\n");
            
            computeNumberOfZeroes(nzeroes, error, x, false);
            result.append("there are approximately ");
            appendMPFRvalue(this.result, nzeroes, 3);
            result.append(" zeroes (with error ");
            appendMPFRvalue(this.result, error, 3);
            result.append(") less than " + this.upperBound + "\n");
            
            MPFR.mpfr_set_si(x, this.lowerBound, mpfr_rnd_t.MPFR_RNDN);
            computeNumberOfZeroes(nzeroes, error, x, true);
            result.append("there are approximately ");
            appendMPFRvalue(this.result, nzeroes, 3);
            result.append(" zeroes (with error ");
            appendMPFRvalue(this.result, error, 3);
            result.append(") less than " + this.lowerBound + "\n");

            publishProgress(-1);
            zeroStr = new StringBuffer[numGramPoints + 1];
            zeroStr[0] = new StringBuffer();
            MPFR.mpfr_set_si(x, this.lowerBound, mpfr_rnd_t.MPFR_RNDN);
    	    if (computeZetaZero(0, x, gramPoints[0], nterms)) {
    		publishProgress(0, 1);
    	    } else {
    		publishProgress(0, 0);
    	    }
            for (i = 1; i < numGramPoints; i++) {
        	zeroStr[i] = new StringBuffer();
        	if (computeZetaZero(i, gramPoints[i - 1], gramPoints[i], nterms)) {
        	    publishProgress(i, 1);
        	} else {
        	    publishProgress(i, 0);
        	}
                if (isCancelled()) {
                    throw new Exception("Task cancelled");
                }
            }
            zeroStr[numGramPoints] = new StringBuffer();
            MPFR.mpfr_set_si(x, this.upperBound, mpfr_rnd_t.MPFR_RNDN);
    	    if (computeZetaZero(numGramPoints, gramPoints[numGramPoints - 1], x, nterms)) {
    		publishProgress(numGramPoints, 1);
    	    } else {
    		publishProgress(numGramPoints, 0);
    	    }
        }
	catch (MPFRException e) {
	    result.append("MPFRException[" + e.getCode() + "] " + e.getMessage());
	    rc = -1;
	}
	catch (GMPException e) {
	    result.append("GMPException[" + e.getCode() + "] " + e.getMessage());
	    rc = -1;
	}
	catch (Exception e) {
	    result.append(e.toString());
	    rc = -1;
	}
        return Integer.valueOf(rc);
    }
    
    protected void onPostExecute(Integer result)
    {
	try {
	    MPFR.mpfr_free_cache();
	}
	catch (MPFRException e) {
	    uinterface.display(TAG + ": " + e.getMessage());
	}
	uinterface.display("number of zeroes found: " + this.numZeroes);
	Log.d(TAG, "number of zeroes found: " + this.numZeroes);
	if (this.result.length() > 0) {
	    uinterface.display(this.result.toString());
	    Log.d(TAG, this.result.toString());
	}
	Log.d(TAG, "number of theta evaluations: " + thetaEvals);
	Log.d(TAG, "number of zeta evaluations: " + zetaEvals);
    }
    
    protected void onPreExecute()
    {
        uinterface.display(TAG);
        uinterface.display("MPFR version " + MPFR.getVersion());
        Log.d(TAG, "MPFR version " + MPFR.getVersion());
        uinterface.display("precision is " + Integer.toString(this.precision) + " bits");
        Log.d(TAG, "precision is " + Integer.toString(this.precision) + " bits");
    }

    protected void onProgressUpdate(Integer... progress)
    {
	int i = progress[0];
	int j;
	String str;
	if (i == -1) {
	    uinterface.display(this.result.toString());
	    Log.d(TAG, this.result.toString());
	    this.result.setLength(0);
	} else {
	    j = progress[1];
	    if (j != 0) {
	        str = "zero[" + Integer.toString(i) + "]= 0.5 + (" + this.zeroStr[i].toString() + ")i";
	    } else {
	        str = "no zero found in " + this.zeroStr[i].toString();
	    }
	    uinterface.display(str);
	    Log.d(TAG, str);
	}
    }

}