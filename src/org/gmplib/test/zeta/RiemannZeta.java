package org.gmplib.test.zeta;

import org.gmplib.gmpjni.MPFR;
import org.gmplib.gmpjni.MPFR.mpfr_t;
import org.gmplib.gmpjni.MPFR.mpfr_rnd_t;
import org.gmplib.gmpjni.MPFR.MPFRException;

public interface RiemannZeta {
    
    void setPrecision(int precision) throws MPFRException, Exception;
    void initBasicFunctions(int n, mpfr_t tmin, mpfr_t tmax) throws MPFRException;
    void initCoefficients(int n) throws MPFRException, Exception;
    void initThetaCoefficients() throws MPFRException;
    void initGramPoints(int lb, int ub, int nt) throws MPFRException;
    int getNumberOfGramPoints();
    int getGramBaseIndex();
    mpfr_t getGramPoint(int i) throws MPFRException;
    void Z(mpfr_t r, mpfr_t t, int n, int nt) throws MPFRException, Exception;
    boolean computeZetaZero(int i, mpfr_t a0, mpfr_t b0, int n, int nt) throws MPFRException, Exception;
    int getNumberOfZeroes();
    int computeNumberOfTerms(mpfr_t tmin, mpfr_t tmax) throws MPFRException, Exception;
    int computeNumberOfThetaTerms(mpfr_t tmin) throws MPFRException;
    long computeNumberOfZeroes(mpfr_t error, mpfr_t t) throws MPFRException, Exception;
    String getZeroStr(int i);
    int getNumberOfZEvals();
    int getNumberOfThetaEvals();

}
