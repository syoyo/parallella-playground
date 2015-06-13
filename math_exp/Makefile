ESDK=${EPIPHANY_HOME}
ELIBS=-L${ESDK}/tools/host/lib
EINCS=-I${ESDK}/tools/host/include
ELDF=${ESDK}/bsps/current/fast.ldf
CROSS_PREFIX=

# How much the host do usleep() to wait a result from e-core?
# Larger value -> longer test time, but can compute much accurate relative error.
WAIT_MICROSECONDS=100000

all:
	echo Build HOST side application
	${CROSS_PREFIX}gcc host.c -o test -DWAIT_MICROSECONDS=${WAIT_MICROSECONDS} ${EINCS} ${ELIBS} -le-hal -lm -le-loader -lpthread
	e-gcc -O3 -g -T ${ELDF} -std=c99 -DFMATH_EXP_TEST=1 -DWAIT_MICROSECONDS=${WAIT_MICROSECONDS} e_fast_exp.c -o e_fast_exp_test.elf -fsingle-precision-constant -mno-soft-cmpsf -mcmove -mfp-mode=truncate -le-lib -lm -ffast-math
	e-objcopy --srec-forceS3 --output-target srec e_fast_exp_test.elf e_fast_exp_test.srec

.PHONY: test
