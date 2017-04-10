/*****************************************************************************
 *   Copyright 2016, 2017 Andy Quick
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *****************************************************************************/
package org.gmplib.test.zeta;

import android.content.Context;
import android.os.AsyncTask;
import android.util.Log;

import org.gmplib.gmpjni.GMP;
import org.gmplib.gmpjni.GMP.GMPException;
import org.gmplib.gmpjni.MPFR;
import org.gmplib.gmpjni.MPFR.mpfr_t;
import org.gmplib.gmpjni.MPFR.mpfr_rnd_t;
import org.gmplib.gmpjni.MPFR.MPFRException;

/**
 * Asynchronous task to compute zeroes of the Riemann zeta function on the critical line.
 * Zeroes are computed on the line 0.5 + i*t for t in [lowerBound, upperBound] by
 * looking for sign changes of Z(t) = exp(i*theta(t)) * zeta(0.5 + i*t) between Gram
 * points in the interval. 
 *
 */
public class Zeta_Task extends AsyncTask<Integer, Integer, Integer>
{

    private static final String TAG = "Zeta_Task";
    private UI uinterface;
    private Context context;
    private StringBuffer result;
    private int precision;
    private int digits; // precision in decimal digits
    private int lowerBound;
    private int upperBound;
    private long elapsedTime;

    private RiemannZeta zetaComputer;
    
    /**
     * Constructor: initialize member variables, get precision from UI, initialize constants.
     */
    public Zeta_Task(UI ui, Context ctx)
        throws MPFRException, Exception
    {
	this.uinterface = ui;
	this.context = ctx;
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
	this.lowerBound = 0;
	this.upperBound = 0;
	this.zetaComputer = null;
    }

    private static int min(int x, int y)
    {
	return (x < y ? x : y);
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
     * Do the work of computing the zeroes.
     */
    protected Integer doInBackground(Integer... params)
    {
	int i;
        int rc = -1;
        int nterms = 20;
        int nthetaterms = 3;
        long t1;
        long t2;
        
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
            if (this.lowerBound >= 200) {
        	result.append("Using Riemann-Siegel formula\n");
                zetaComputer = new RiemannSiegelZeta();
            } else if (this.upperBound <= 200) {
        	result.append("Using Riemann zeta basic formula\n");
                zetaComputer = new RiemannZetaBasic();
            } else {
        	throw new Exception("lower and upper bounds must both be less than 200 or greater than 200");
            }
            zetaComputer.setPrecision(precision);
            mpfr_t a = new mpfr_t();
            mpfr_t b = new mpfr_t();
            long nzeroes;
            mpfr_t error = new mpfr_t();
            this.digits = (int)((double)(precision - 2)*0.301029995664 /* log10(2) */);

            zetaComputer.initThetaCoefficients();
            MPFR.mpfr_set_si(a, this.lowerBound, mpfr_rnd_t.MPFR_RNDN);
            MPFR.mpfr_set_si(b, this.upperBound, mpfr_rnd_t.MPFR_RNDN);
            nthetaterms = zetaComputer.computeNumberOfThetaTerms(a);
            zetaComputer.initGramPoints(this.lowerBound, this.upperBound, min(3, nthetaterms));
            result.append("Gram Points\n");
            for (i = 0; i < zetaComputer.getNumberOfGramPoints(); i++) {
        	result.append("    g(" + Integer.toString(zetaComputer.getGramBaseIndex() + i) + ")=");
        	appendMPFRvalue(this.result, zetaComputer.getGramPoint(i), 5);
        	result.append("\n");
            }
            nterms = zetaComputer.computeNumberOfTerms(a, b);
            zetaComputer.initBasicFunctions(nterms, a, b);
            zetaComputer.initCoefficients(nterms);
            if (zetaComputer instanceof RiemannSiegelZeta) {
                result.append("using " + nterms + " low-order terms in Z sum\n");
            } else {
                result.append("using " + nterms + " terms in zeta sums\n");        	
            }
            result.append("using " + nthetaterms + " low-order terms in theta\n");
            
            nzeroes = zetaComputer.computeNumberOfZeroes(error, b);
            if (MPFR.mpfr_cmp_si(error, 0) == 0) {
                result.append("there are ");
                result.append(nzeroes);
                result.append(" zeroes");        	
            } else {
                result.append("there are approximately ");
                result.append(nzeroes);
                result.append(" zeroes (with error ");
                appendMPFRvalue(this.result, error, 3);
                result.append(")");
            }
            result.append(" less than " + this.upperBound + "\n");
            
            nzeroes = zetaComputer.computeNumberOfZeroes(error, a);
            if (MPFR.mpfr_cmp_si(error, 0) == 0) {
                result.append("there are ");
                result.append(nzeroes);
                result.append(" zeroes");        	
            } else {
                result.append("there are approximately ");
                result.append(nzeroes);
                result.append(" zeroes (with error ");
                appendMPFRvalue(this.result, error, 3);
                result.append(")");
            }
            result.append(" less than " + this.lowerBound + "\n");

            publishProgress(-1);
            t1 = System.currentTimeMillis();
    	    if (zetaComputer.computeZetaZero(0, a, zetaComputer.getGramPoint(0), nterms, nthetaterms)) {
    		publishProgress(0, 1);
    	    } else {
    		publishProgress(0, 0);
    	    }
            int n = zetaComputer.getNumberOfGramPoints();
            for (i = 1; i < n; i++) {
        	if (zetaComputer.computeZetaZero(i, zetaComputer.getGramPoint(i - 1), zetaComputer.getGramPoint(i), nterms, nthetaterms)) {
        	    publishProgress(i, 1);
        	} else {
        	    publishProgress(i, 0);
        	}
                if (isCancelled()) {
                    throw new Exception("Task cancelled");
                }
            }
    	    if (zetaComputer.computeZetaZero(n, zetaComputer.getGramPoint(n - 1), b, nterms, nthetaterms)) {
    		publishProgress(n, 1);
    	    } else {
    		publishProgress(n, 0);
    	    }
    	    t2 = System.currentTimeMillis();
    	    this.elapsedTime = t2 - t1;
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
	    /***
	    StackTraceElement[] st = e.getStackTrace();
	    for (int m = 0; m < st.length; m++) {
		result.append("\n");
		result.append(st[m].toString());
	    }
	    ***/
	    rc = -1;
	}
        return Integer.valueOf(rc);
    }
    
    /**
     * Post-execution work.
     */
    protected void onPostExecute(Integer result)
    {
	try {
	    MPFR.mpfr_free_cache();
	}
	catch (MPFRException e) {
	    uinterface.display(TAG + ": " + e.getMessage());
	}
	if (this.zetaComputer != null) {
	    uinterface.display("number of zeroes found: " + this.zetaComputer.getNumberOfZeroes());
	    Log.d(TAG, "number of zeroes found: " +  this.zetaComputer.getNumberOfZeroes());
	}
	if (this.result.length() > 0) {
	    uinterface.display(this.result.toString());
	    Log.d(TAG, this.result.toString());
	}
	if (this.zetaComputer != null) {
	    Log.d(TAG, "number of theta evaluations: " + this.zetaComputer.getNumberOfThetaEvals());
	    Log.d(TAG, "number of Z evaluations: " + this.zetaComputer.getNumberOfZEvals());
	    Log.d(TAG, "elapsed time: " + this.elapsedTime + " milliseconds");
	}
    }
    
    /**
     * Pre-execution work.
     */
    protected void onPreExecute()
    {
	String str;
        uinterface.display(TAG);
        uinterface.display(context.getString(R.string.gram_note));
        str = "MPFR " + context.getString(R.string.version) + " " + MPFR.getVersion();
        uinterface.display(str);
        Log.d(TAG, str);
        str = context.getString(R.string.precision_is) + " " + Integer.toString(this.precision) + " " + context.getString(R.string.bits);
        uinterface.display(str);
        Log.d(TAG, str);
    }

    /**
     * Update the UI with each zero as it is found.
     */
    protected void onProgressUpdate(Integer... progress)
    {
	int i = progress[0];
	int j;
	String str;
	StringBuffer sb = new StringBuffer();
	if (i == -1) {
	    uinterface.display(this.result.toString());
	    Log.d(TAG, this.result.toString());
	    this.result.setLength(0);
	} else {
	    sb.append("[");
	    sb.append(Integer.toString(i));
	    sb.append("] ");
	    j = progress[1];
	    if (j != 0) {
	        sb.append("0.5 + (");
	        sb.append(this.zetaComputer.getZeroStr(i));
	        sb.append(")i");
	    } else {
	        sb.append("no zero found in ");
	        sb.append(this.zetaComputer.getZeroStr(i));
	    }
	    str = sb.toString();
	    uinterface.display(str);
	    Log.d(TAG, str);
	}
    }

}
