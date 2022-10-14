//
// Created by mlevorato on 09/06/18.
//

#ifndef FLOWSHOP_SOLVER_NUMERICUTIL_H
#define FLOWSHOP_SOLVER_NUMERICUTIL_H

#include <cmath>
#include <cassert>

// number precision macros obtained from Test project (http://scip.zib.de/doc/html/def_8h.php)
#define 	REALABS(x)   (fabs(x))

#define 	EPSEQ(x, y, eps)   (REALABS((x)-(y)) <= (eps))

#define 	EPSLT(x, y, eps)   ((x)-(y) < -(eps))

#define 	EPSLE(x, y, eps)   ((x)-(y) <= (eps))

#define 	EPSGT(x, y, eps)   ((x)-(y) > (eps))

#define 	EPSGE(x, y, eps)   ((x)-(y) >= -(eps))

#define 	EPSZ(x, eps)   (REALABS(x) <= (eps))

#define 	EPSP(x, eps)   ((x) > (eps))

#define 	EPSN(x, eps)   ((x) < -(eps))

#define 	EPSFLOOR(x, eps)   (floor((x)+(eps)))

#define 	EPSCEIL(x, eps)   (ceil((x)-(eps)))

#define 	EPSROUND(x, eps)   (ceil((x)-0.5+(eps)))

#define 	EPSFRAC(x, eps)   ((x)-EPSFLOOR(x,eps))

#define 	EPSISINT(x, eps)   (EPSFRAC(x,eps) <= (eps))

#define 	Test_DEFAULT_INFINITY   1e+20

#define 	Test_DEFAULT_EPSILON   1e-09

#define 	Test_DEFAULT_SUMEPSILON   1e-06

#define 	Test_DEFAULT_HUGEVAL   1e+15

#define 	Test_MAXEPSILON   1e-03

#define 	Test_MINEPSILON   1e-20


class NumericUtil {
public:

    static double num_epsilon; /**< absolute values smaller than this are considered zero */
    static double num_infinity;
    static double num_hugeval;

    /** checks, if value is in range epsilon of 0.0 */
   static bool TestsetIsZero(
      double             val                 /**< value to process */
      )
   {
      return EPSZ(val, num_epsilon);
   }

    /** checks, if value is greater than epsilon */
    static bool TestsetIsPositive(
            double             val                 /**< value to process */
    )
    {
       return EPSP(val, num_epsilon);
    }

    /** checks, if value is lower than -epsilon */
    static bool TestsetIsNegative(
            double             val                 /**< value to process */
    )
    {
       return EPSN(val, num_epsilon);
    }

    /** checks, if value is (positive) infinite */
    static bool TestsetIsInfinity(
            
            double             val                 /**< value to be compared against infinity */
    )
    {
        return (val >= num_infinity);
    }

    /** checks, if value is huge and should be handled separately (e.g., in activity computation) */
    static bool TestsetIsHugeValue(
            
            double             val                 /**< value to be checked whether it is huge */
    )
    {
        return (val >= num_hugeval);
    }

    /** checks, if values are in range of epsilon */
    static bool TestsetIsEQ(
            
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        return EPSEQ(val1, val2, num_epsilon);
    }

    /** checks, if val1 is (more than epsilon) lower than val2 */
    static bool TestsetIsLT(
            
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        return EPSLT(val1, val2, num_epsilon);
    }

    /** checks, if val1 is not (more than epsilon) greater than val2 */
    static bool TestsetIsLE(
            
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        return EPSLE(val1, val2, num_epsilon);
    }

    /** checks, if val1 is (more than epsilon) greater than val2 */
    static bool TestsetIsGT(
            
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        return EPSGT(val1, val2, num_epsilon);
    }

    /** checks, if val1 is not (more than epsilon) lower than val2 */
    static bool TestsetIsGE(
            
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        return EPSGE(val1, val2, num_epsilon);
    }

    /** checks, if val1 is (more than epsilon) lower than val2 */
    static bool TestisLT(
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        /* avoid to compare two different infinities; the reason for that is
         * that such a comparison can lead to unexpected results */
        assert( ((!TestisInfinity(val1) || !TestisInfinity(val2))
                 && (!TestisInfinity(-val1) || !TestisInfinity(-val2)))
                || val1 == val2 );    /*lint !e777*/

        return TestsetIsLT(val1, val2);
    }

/** checks, if val1 is not (more than epsilon) greater than val2 */
    static bool TestisLE(
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        /* avoid to compare two different infinities; the reason for that is
         * that such a comparison can lead to unexpected results */
        assert( ((!TestisInfinity(val1) || !TestisInfinity(val2))
                 && (!TestisInfinity(-val1) || !TestisInfinity(-val2)))
                || val1 == val2 );    /*lint !e777*/

        return TestsetIsLE(val1, val2);
    }

/** checks, if val1 is (more than epsilon) greater than val2 */
    static bool TestisGT(
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        /* avoid to compare two different infinities; the reason for that is
         * that such a comparison can lead to unexpected results */
        assert( ((!TestisInfinity(val1) || !TestisInfinity(val2))
                 && (!TestisInfinity(-val1) || !TestisInfinity(-val2)))
                || val1 == val2 );    /*lint !e777*/

        return TestsetIsGT(val1, val2);
    }

/** checks, if val1 is not (more than epsilon) lower than val2 */
    static bool TestisGE(
            double             val1,               /**< first value to be compared */
            double             val2                /**< second value to be compared */
    )
    {
        /* avoid to compare two different infinities; the reason for that is
         * that such a comparison can lead to unexpected results */
        assert( ((!TestisInfinity(val1) || !TestisInfinity(val2))
                 && (!TestisInfinity(-val1) || !TestisInfinity(-val2)))
                || val1 == val2 );    /*lint !e777*/

        return TestsetIsGE(val1, val2);
    }

    /** returns value treated as infinity */
    static double Testinfinity(
    )
    {
        return 1e+20;
    }

/** checks, if value is (positive) infinite */
    static bool TestisInfinity(
            double             val                 /**< value to be compared against infinity */
    )
    {
        return TestsetIsInfinity(val);
    }

/** checks, if value is huge and should be handled separately (e.g., in activity computation) */
    static bool TestisHugeValue(
            double             val                 /**< value to be checked whether it is huge */
    )
    {
        return TestsetIsHugeValue(val);
    }

/** returns the minimum value that is regarded as huge and should be handled separately (e.g., in activity
 *  computation)
 */
    static double TestgetHugeValue(
    )
    {
        return 1e+15;
    }

/** checks, if value is in range epsilon of 0.0 */
    static bool TestisZero(
            double             val                 /**< value to process */
    )
    {
        return TestsetIsZero(val);
    }

/** checks, if value is greater than epsilon */
    static bool TestisPositive(
            double             val                 /**< value to process */
    )
    {
        return TestsetIsPositive(val);
    }

/** checks, if value is lower than -epsilon */
    static bool TestisNegative(
            double             val                 /**< value to process */
    )
    {
        return TestsetIsNegative(val);
    }


};


#endif //FLOWSHOP_SOLVER_NUMERICUTIL_H
