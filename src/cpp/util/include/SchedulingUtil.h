//
// Created by Mario Costa Levorato Junior on 2019-01-11.
// Source: http://new.zsd.iiar.pwr.wroc.pl/files/sources/sax.cpp
//

#ifndef FLOWSHOP_SOLVER_SCHEDULINGUTIL_H
#define FLOWSHOP_SOLVER_SCHEDULINGUTIL_H

#include <cstdlib>
#include <cstdio>
// #include <conio.h>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include "../../problem/GlobalTypes.h"

// FIXME adjust the maximum number of machines and jobs according to larger problem sizes!
#define NJOBS    1002
#define MACHINES    502

namespace util {

    using namespace std;

    class TwoMachineLB {
    public:
        TwoMachineLB(std::pair<int, int> mach_pair, double lb) : machine_pair(mach_pair), value(lb) {
        }

        std::pair<int, int> machine_pair;
        double value;
    };

    class SchedulingUtil {
    public:

        template<class number>
        int sign(number x) { return (x <= 0) ? -1 : 1; }

        // swaps two elements
        template<class vect>
        static void swap(vect *a, vect *b) {
            vect c = *a;
            *a = *b;
            *b = c;
        }

        /**
         * quicksort; returns permutation such that a[pi[i]]<=a[pi[i+1]] for a[n..m]
         *
         * @tparam vect
         * @param n
         * @param m
         * @param a
         * @param pi
         */
        template<class vect>
        static void sort(int n, int m, vect a[], int pi[]) {
            int i, j;
            vect x;
            if (m <= n) return;
            i = n;
            j = m;
            x = a[pi[(i + j) / 2]];
            do {
                while (a[pi[i]] < x) i++;
                while (x < a[pi[j]]) j--;
                if(i >= j) break;
                if (i < j) {  // correcting previous bug
                    swap(pi + i, pi + j);
                    i++;
                    j--;
                }
            } while (i < j);
            sort(n, i-1, a, pi);  // correcting previous bug
            sort(j+1, m, a, pi);
        }

        /**
         * quicksort; returns permutation such that a[pi[i]]>=a[pi[i+1]] for a[n..m]
         *
         * @tparam vect
         * @param n
         * @param m
         * @param a
         * @param pi
         */
        template<class vect>
        static void _sort(int left, int right, vect a[], int pi[]) {
            if(left >= right) return;
            int i = left, j = right;
            vect pivot = a[pi[i]];
            int tmp = 0;
            do {
                while(a[pi[i]] > pivot) i++;
                while(pivot > a[pi[j]]) j--;
                if(i >= j) break;
                tmp = pi[i]; pi[i] = pi[j]; pi[j] = tmp;
                i++; j--;
            } while (i < j);
            _sort(left, i-1, a, pi);  // correcting previous bug
            _sort(j+1, right, a, pi);
        }

        // --- nondecreasing sorting
        static void quicksort(long a[], int pi[], int n, int m)
        {
            int i,j,p; long x;

            if (m <= n) return;
            i=n; j=m; x=a[pi[(i+j)>>1]];
            do
            {  while (a[pi[i]] < x) i++;
                while (x < a[pi[j]]) j--;
                if (i <= j) { p=pi[j]; pi[j]=pi[i]; pi[i]=p; i++; j--; }
            }  while (i < j);
            quicksort(a,pi,n,j); quicksort(a,pi,i,m);
        }

        /**
         * Cmax for permutation pi
         * @param n
         * @param m
         * @param p
         * @param pi
         * @return
         */
        static long Cmax(int n, int m, long p[][NJOBS + 1], int pi[]) {
            int i, j, k;
            long C[MACHINES + 1][NJOBS + 1];

            k = pi[1];
            C[1][1] = p[1][k];
            for (i = 2; i <= m; i++) C[i][1] = C[i - 1][1] + p[i][k];

            for (j = 2; j <= n; j++) {
                k = pi[j];
                C[1][j] = C[1][j - 1] + p[1][k];
                for (i = 2; i <= m; i++) C[i][j] = std::max(C[i][j - 1], C[i - 1][j]) + p[i][k];
            }
            return C[m][n];
        }

        /**
         * Cmax for permutation pi, using ublas::matrix as input for proc times
         *
         * @param n
         * @param m
         * @param p
         * @param pi
         * @return
         */
        static double Cmax(int n, int m, const ublas::matrix<double> &p, int pi[]) {
            int i, j, k;
            double C[MACHINES + 1][NJOBS + 1];

            k = pi[1];
            C[1][1] = p(1, k);
            for (i = 2; i <= m; i++) C[i][1] = C[i - 1][1] + p(i, k);

            for (j = 2; j <= n; j++) {
                k = pi[j];
                C[1][j] = C[1][j - 1] + p(1, k);
                for (i = 2; i <= m; i++) C[i][j] = std::max(C[i][j - 1], C[i - 1][j]) + p(i, k);
            }
            return C[m][n];
        }

        static long power(int x, int n) {
            int i;
            long s = 1;
            for (i = 1; i <= n; i++) s *= x;
            return s;
        }

        static long eps(int x, int n) {
            int i;
            long s = 0;
            for (i = 1; i <= n; i++) s += power(x, i - 1);
            return s;
        }

    };
}

#endif //FLOWSHOP_SOLVER_SCHEDULINGUTIL_H
