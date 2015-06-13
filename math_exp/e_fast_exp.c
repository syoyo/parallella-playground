#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "e_lib.h"

// GCC
#define RESTRICT __restrict__

extern float expapprox(float val);
extern void expapprox4(float *RESTRICT dst, const float *src);

//
// -------------------------------------------------------------------------------------
// Fast approximate exp() function for Epiphany.
//
// Copyright 2015 Syoyo Fujita <syoyo@lighttransport.com>
// Licensed under Apache 2.0 License.
//
//
// <<Epiphany>>
//   -O3 -std=c99 -fsingle-precision-constant -mno-soft-cmpsf -mcmove
//   -mfp-mode=truncate
//
//                          Enable range check       Disable range check
//  -------------+--------+-------------------------+-------------------------------
//  expapprox()  | scalar | 74 cycles               | 54 cycles
//  expapprox4() | 4 SIMD | 32 cycles(128 in total) | 25 cycles(90 in total)
//
//
//  - Reference
//
//  -------------+--------+---------------------------------------------------------
//  expf()       | scalar | 141633 cycles(tooooo slow because of SW
//  implementation?)
//
//  - Difference
//
//    [-30.0, 30.0] Relative diff: ave = 0.000002, min = 0.000000, max =
//    0.000008
//
//  - Note
//    * -mfp-mode=truncate may loose some precison, but emits more optimal
//    assembly.
//    * -mno-soft-compsf required to not emit software compare function, which
//    breaks
//      some IEEE 754 compliance, but should be totally OK for approximate math
//      functions.
//
//

// If you are sure that the input value for exp() is within [-88.0, 88.0],
// you can disable range check, which results in faster evaluation.
#define FMATH_EXP_DISABLE_RANGE_CHECK (0)

// Based on http://gallium.inria.fr/blog/fast-vectorizable-math-approx/

/* Relative error bounded by 1e-5 for normalized outputs
   Returns invalid outputs for nan inputs */
float expapprox(float val)
{
	/* Workaround a lack of optimization in gcc */
	const float exp_cst1 = 2139095040.f;
	const float exp_cst2 = 0.f;

	union {
		int i;
		float f;
	} xu, xu2;
	float val2, val3, val4, b;
	int val4i;
	val2 = 12102203.1615614f * val + 1065353216.f;
#if FMATH_EXP_DISABLE_RANGE_CHECK
	val4 = val2;
#else
	val3 = val2 < exp_cst1 ? val2 : exp_cst1;
	val4 = val3 > exp_cst2 ? val3 : exp_cst2;
#endif
	val4i = (int)val4;
	xu.i = val4i & 0x7F800000;
	xu2.i = (val4i & 0x7FFFFF) | 0x3F800000;
	b = xu2.f;
	return xu.f *
	       (0.509964287281036376953125f +
		b * (0.3120158612728118896484375f +
		     b * (0.1666135489940643310546875f +
			  b * (-2.12528370320796966552734375e-3f +
			       b * 1.3534179888665676116943359375e-2f))));
}

void expapprox4(float *RESTRICT dst, const float *RESTRICT src)
{
	// Manual code expansion of exparrpox() x 4.

	/* Workaround a lack of optimization in gcc */
	const float exp_cst1 = 2139095040.f;
	const float exp_cst2 = 0.f;

	const float kCoeff[5] = {
	    0.509964287281036376953125f, 0.3120158612728118896484375f,
	    0.1666135489940643310546875f, -2.12528370320796966552734375e-3f,
	    1.3534179888665676116943359375e-2f};

	union {
		int i;
		float f;
	} xu_0, xu_1, xu_2, xu_3, xu2_0, xu2_1, xu2_2, xu2_3;

	float val2_0, val2_1, val2_2, val2_3;
	float val3_0, val3_1, val3_2, val3_3;
	float val4_0, val4_1, val4_2, val4_3;
	float b0, b1, b2, b3;
	int val4i_0, val4i_1, val4i_2, val4i_3;

	val2_0 = 12102203.1615614f * src[0] + 1065353216.f;
	val2_1 = 12102203.1615614f * src[1] + 1065353216.f;
	val2_2 = 12102203.1615614f * src[2] + 1065353216.f;
	val2_3 = 12102203.1615614f * src[3] + 1065353216.f;

#if FMATH_EXP_DISABLE_RANGE_CHECK
	val4_0 = val2_0;
	val4_1 = val2_1;
	val4_2 = val2_2;
	val4_3 = val2_3;
#else
	val3_0 = val2_0 < exp_cst1 ? val2_0 : exp_cst1;
	val3_1 = val2_1 < exp_cst1 ? val2_1 : exp_cst1;
	val3_2 = val2_2 < exp_cst1 ? val2_2 : exp_cst1;
	val3_3 = val2_3 < exp_cst1 ? val2_3 : exp_cst1;

	val4_0 = val3_0 > exp_cst2 ? val3_0 : exp_cst2;
	val4_1 = val3_1 > exp_cst2 ? val3_1 : exp_cst2;
	val4_2 = val3_2 > exp_cst2 ? val3_2 : exp_cst2;
	val4_3 = val3_3 > exp_cst2 ? val3_3 : exp_cst2;
#endif

	val4i_0 = (int)val4_0;
	val4i_1 = (int)val4_1;
	val4i_2 = (int)val4_2;
	val4i_3 = (int)val4_3;

	xu_0.i = val4i_0 & 0x7F800000;
	xu_1.i = val4i_1 & 0x7F800000;
	xu_2.i = val4i_2 & 0x7F800000;
	xu_3.i = val4i_3 & 0x7F800000;

	xu2_0.i = (val4i_0 & 0x7FFFFF) | 0x3F800000;
	xu2_1.i = (val4i_1 & 0x7FFFFF) | 0x3F800000;
	xu2_2.i = (val4i_2 & 0x7FFFFF) | 0x3F800000;
	xu2_3.i = (val4i_3 & 0x7FFFFF) | 0x3F800000;

	b0 = xu2_0.f;
	b1 = xu2_1.f;
	b2 = xu2_2.f;
	b3 = xu2_3.f;

	const float c0 = kCoeff[0];
	const float c1 = kCoeff[1];
	const float c2 = kCoeff[2];
	const float c3 = kCoeff[3];
	const float c4 = kCoeff[4];

	dst[0] = xu_0.f * (c0 + b0 * (c1 + b0 * (c2 + b0 * (c3 + b0 * c4))));
	dst[1] = xu_1.f * (c0 + b1 * (c1 + b1 * (c2 + b1 * (c3 + b1 * c4))));
	dst[2] = xu_2.f * (c0 + b2 * (c1 + b2 * (c2 + b2 * (c3 + b2 * c4))));
	dst[3] = xu_3.f * (c0 + b3 * (c1 + b3 * (c2 + b3 * (c3 + b3 * c4))));
}

//
// -------------------------------------------------------------------------------------
//

#if FMATH_EXP_TEST

#ifndef WAIT_MICROSECONDS
#define WAIT_MICROSECONDS (10000)
#endif

#include <stdio.h>

void validateExp(float retDiff[3], float beginValue, float endValue)
{
	int n = WAIT_MICROSECONDS / 250; // 250 = emprically found value.
	union {
		int i;
		float f;
	} bv, ev, it;
	bv.f = beginValue;
	ev.f = endValue;
	float step = (endValue - beginValue) / n;

	it.i = bv.i;
	int count = 0;
	volatile float minDiff = 0.0f;
	volatile float maxDiff = 0.0f;
	volatile float aveDiff = 0.0f;
	// for (; it.i < ev.i; it.i++) { // <-- e-gcc can't compile this loop
	// ...
	float f = beginValue;
	for (f = beginValue; f < endValue; f += step) {

		float ref = expf(f);
		float ret = expapprox(f);
		float diff = fabsf(ref - ret) / ref;

		if (count == 0) {
			minDiff = diff;
			maxDiff = diff;
		} else {
			minDiff = (minDiff > diff) ? diff : minDiff;
			maxDiff = (maxDiff < diff) ? diff : maxDiff;
		}
		aveDiff += diff;
		count++;
	}

	aveDiff /= (float)count;

	retDiff[0] = aveDiff;
	retDiff[1] = minDiff;
	retDiff[2] = maxDiff;
}

void validateExp4(float retDiff[3], float beginValue, float endValue)
{
	int n = WAIT_MICROSECONDS / 250; // 250 = emprically found value.
	union {
		int i;
		float f;
	} bv, ev, it;
	bv.f = beginValue;
	ev.f = endValue;
	float step = (endValue - beginValue) / n;

	it.i = bv.i;
	int count = 0;
	volatile float minDiff = 0.0f;
	volatile float maxDiff = 0.0f;
	volatile float aveDiff = 0.0f;
	// for (; it.i < ev.i; it.i++) { // <-- e-gcc can't compile this loop
	// ...
	float f = beginValue;
	for (f = beginValue; f < endValue; f += 4.0f * step) {

		float ref[4];
		ref[0] = expf(f);
		ref[1] = expf(f + step);
		ref[2] = expf(f + 2 * step);
		ref[3] = expf(f + 3 * step);
		float ret[4];
		float src[4] = {f, f + step, f + 2 * step, f + 3 * step};
		expapprox4(ret, src);

		for (int k = 0; k < 4; k++) {
			float diff = fabs(ref[k] - ret[k]) / (float)ref[k];

			if (count == 0) {
				minDiff = diff;
				maxDiff = diff;
			} else {
				minDiff = (minDiff > diff) ? diff : minDiff;
				maxDiff = (maxDiff < diff) ? diff : maxDiff;
			}
			aveDiff += diff;
			count++;
		}
	}

	aveDiff /= (float)count;

	retDiff[0] = aveDiff;
	retDiff[1] = minDiff;
	retDiff[2] = maxDiff;
}

char outbuf[4096] SECTION("shared_dram");
int main(void)
{
	e_coreid_t coreid;
	unsigned int i;
	unsigned int num;
	unsigned int time_p;
	unsigned int time_c;
	unsigned int time_compare;
	unsigned *mailbox;
	unsigned int temp;

	float volatile in_sin;
	float volatile in_cos;
	float volatile in_sqt;
	float volatile in_ceil;
	float volatile in_log;
	float volatile in_exp;
	float volatile in_exp1;
	float volatile in_exp2;
	float volatile in_exp3;
	float volatile in_exp4;
	float volatile in_exp5;
	float volatile in_exp6;
	float volatile in_exp7;
	float re_f0, re_f1, re_f2, re_f3, re_f4, re_f5;

	in_exp = 1.88f;
	in_exp1 = 2.88f;
	in_exp2 = 3.88f;
	in_exp3 = 4.88f;
	in_exp4 = 4.88f;
	in_exp5 = 5.88f;
	in_exp6 = 6.88f;
	in_exp7 = 7.88f;
	mailbox = (unsigned *)0x6000;
	mailbox[0] = 0;
	mailbox[1] = 0xFFFFFFFF;
	mailbox[2] = 0xFFFFFFFF;
	mailbox[3] = 0xFFFFFFFF;

	// Who am I? Query the CoreID from hardware.
	coreid = e_get_coreid();
	sprintf(outbuf, "");

	const float in_exp_arr[8] = {in_exp,  in_exp1, in_exp2, in_exp3,
				     in_exp4, in_exp5, in_exp6, in_exp7};
	float out_exp_arr[8];

	// Get time waste on functions
	e_ctimer_set(E_CTIMER_0, E_CTIMER_MAX);
	time_p = e_ctimer_start(E_CTIMER_0, E_CTIMER_CLK);
	time_c = e_ctimer_get(E_CTIMER_0);
	e_ctimer_stop(E_CTIMER_0);
	time_compare = time_p - time_c;

	// expf() reference
	volatile unsigned int exp_ref_clocks = 0;
	{
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		volatile float ret = expf(in_exp);

		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		exp_ref_clocks = time_p - time_c - time_compare;

		sprintf(outbuf + strlen(outbuf), "\nThe clock cycle count for "
						 "\"expf()\" (reference) is "
						 "%d.\n",
			exp_ref_clocks);
	}

	// expapprox
	if (1) {
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		// re_f4 = expapprox(in_exp);
		// re_f4 = fmath_exp(in_exp);
		// fmath_exp4(out_exp_arr, in_exp_arr);
		// fmath_exp8(out_exp_arr, in_exp_arr);
		// volatile float ret[4];
		volatile float ret = expapprox(in_exp);
		// ret[1] = expapprox(in_exp+0.2);
		// ret[2] = expapprox(in_exp+0.3);
		// ret[3] = expapprox(in_exp+0.4);
		// expapprox(out_exp_arr, in_exp_arr);
		// expapprox8(out_exp_arr, in_exp_arr);
		// re_f4 = expapprox(in_exp) + expapprox(in_exp1) +
		// expapprox(in_exp2) +
		// expapprox(in_exp3); //ret[0] + ret[1] + ret[2] + ret[3];
		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		// re_f4 = out_exp_arr[0] + out_exp_arr[1] + out_exp_arr[2] +
		// out_exp_arr[3]
		// + out_exp_arr[4] + out_exp_arr[5] + out_exp_arr[6] +
		// out_exp_arr[7];

		temp = time_p - time_c - time_compare;

		// volatile float ref = expf(in_exp);
		// volatile float my = expapprox(in_exp);
		// volatile float diff = fabsf(ref - my);
		sprintf(outbuf + strlen(outbuf),
			"\nThe clock cycle count for \"expapprox()\" is %d.\n",
			temp);
		// sprintf(outbuf+strlen(outbuf)  , "\ndiff: %d.\n", (int)diff);
	}

	// expapprox4
	if (1) {
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		expapprox4(out_exp_arr, in_exp_arr);
		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		// prevent compiler dead code optimization
		volatile float ret = out_exp_arr[0] + out_exp_arr[1] +
				     out_exp_arr[2] + out_exp_arr[3];
		// re_f4 = out_exp_arr[0] + out_exp_arr[1] + out_exp_arr[2] +
		// out_exp_arr[3]
		// + out_exp_arr[4] + out_exp_arr[5] + out_exp_arr[6] +
		// out_exp_arr[7];

		temp = time_p - time_c - time_compare;

		// volatile float ref = expf(in_exp);
		// volatile float my = expapprox(in_exp);
		// volatile float diff = fabsf(ref - my);
		sprintf(outbuf + strlen(outbuf), "\nThe clock cycle count for "
						 "\"expapprox4()\" is %d (/4 = "
						 "%d).\n",
			temp, temp / 4);
		// sprintf(outbuf+strlen(outbuf)  , "\ndiff: %d.\n", (int)diff);
	}

	// Validation
	{
		float diffs[3];
		// validateExp(diffs, -30.0f, 30.0f);
		validateExp4(diffs, -3.0f, 3.0f);

		mailbox[0] = 1;
		mailbox[1] = *((unsigned int *)&diffs[0]); // ave
		mailbox[2] = *((unsigned int *)&diffs[1]); // min
		mailbox[3] = *((unsigned int *)&diffs[2]); // max
	}

	return EXIT_SUCCESS;
}
#endif
