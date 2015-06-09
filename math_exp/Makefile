all:
	e-gcc -O3 -mno-soft-cmpsf -mfp-mode=truncate -S fmath_exp.c

gen:
	gcc -o tablegen fmath_exp_gen.c -lm
