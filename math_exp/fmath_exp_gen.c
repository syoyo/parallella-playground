//
// TODO: We can reduce table size since its value use only 23 bits for each table element.
//
#include <math.h>
#include <stdio.h>

// base on fmath::exp https://github.com/herumi/fmath/blob/master/fmath.hpp

inline unsigned int mask(int x)
{
	return (1U << x) - 1;
}

typedef union {
	float f;
	unsigned int i;
} fi;

void fmath_exp_gentable(int tableSize) {
	
	const int n = 1 << tableSize;
	int i = 0;

	printf("const unsigned int kFmathExpTable[%d] = {\n  ", 1 << tableSize);

	for (i = 0; i < n; i++) {
		float y = pow(2.0f, (float)i / n);
		fi fi;
		fi.f = y;
		unsigned int value = fi.i & mask(23);

		printf("0x%08x", value);
		if (i != (n-1)) printf(", ");
		if ((i != 0) && (i % 4 == 3)) {
			printf("\n");
			if (i != (n-1)) {
				printf("  ");
			}
		}

	}
	printf("};\n");
}

int main(int argc, char** argv)
{
	int tableSize = 10;
	if (argc > 1) {
		tableSize = atoi(argv[1]);
	}

	if (tableSize > (1 << 14)) {
		tableSize = (1 << 14);
	}

	fmath_exp_gentable(tableSize);
} 

