#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "e_lib.h"
#ifdef __cplusplus
}
#endif


// GCC
#define RESTRICT __restrict__

// rayov = rayorg * rayinvdir
char ray_aabb(float outT[2], float maxT, const float bbox[2][3], const float rayov[3], const float rayinvdir[3], const char raydirsign[3]) {
#if 1
	const float min_x = bbox[raydirsign[0]^1][0];	
	const float min_y = bbox[raydirsign[1]^1][1];	
	const float min_z = bbox[raydirsign[2]^1][2];	

	const float max_x = bbox[raydirsign[0]][0];	
	const float max_y = bbox[raydirsign[1]][1];	
	const float max_z = bbox[raydirsign[2]][2];	

	const float tmin_x = min_x * rayinvdir[0] + rayov[0];
	const float tmax_x = max_x * rayinvdir[0] + rayov[0];

	const float tmin_y = min_y * rayinvdir[1] + rayov[1];
	const float tmax_y = max_y * rayinvdir[1] + rayov[1];

	float tmin = (tmin_x > tmin_y) ? tmin_x : tmin_y;
	float tmax = (tmax_x < tmax_y) ? tmax_x : tmax_y;

	const float tmin_z = min_z * rayinvdir[2] + rayov[2];
	const float tmax_z = max_z * rayinvdir[2] + rayov[2];

	tmin = (tmin > tmin_z) ? tmin : tmin_z;
	tmax = (tmax < tmax_z) ? tmax : tmax_z;

	// Write out tmin/tmax anyway for the performane.
	outT[0] = tmin;		
	outT[1] = tmax;		

	char hit = (tmax > 0.0f) && (tmin <= tmax) && (tmin <= maxT);
#endif

	return 0; // false
}

#if RAYTRACE_TEST

#ifndef WAIT_MICROSECONDS
#define WAIT_MICROSECONDS (10000)
#endif

#define TEST_NUM	  (2)

#include <stdio.h>

char outbuf[4096] SECTION("shared_dram");
int main(int argc, char **argv)
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

	float outT[2];
	volatile float maxT = 10.0f;
	const float bbox[2][3] = {{-1.0f+argc, -1.0f+argc, -1.0f+argc}, {1.0f, 1.0f, 1.0f}};
	const float rayov[3] = {0.0f, 0.0f, 0.0f};
	const float rayinvdir[3] = {1.0f+argc, 2.0f, 3.0f};
	char  raydirsign[3] = {0,0,0};

	volatile unsigned int code_clocks = 0;
	{
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		volatile char hit = ray_aabb(outT, maxT, bbox, rayov, rayinvdir, raydirsign);

		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		code_clocks = time_p - time_c - time_compare;

		sprintf(outbuf + strlen(outbuf), "\nThe clock cycle count for "
						 "\"ray_aabb()\" is "
						 "%d.\n",
			code_clocks);

		mailbox[0] = 1;
	}

	return EXIT_SUCCESS;
}
#endif
