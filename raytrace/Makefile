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
	e-g++ -O3 -g -T ${ELDF} -DRAYTRACE_TEST=1 -DWAIT_MICROSECONDS=${WAIT_MICROSECONDS} e_raytrace.cc -o e_raytrace.elf -fsingle-precision-constant -mno-soft-cmpsf -mcmove -mfp-mode=truncate -le-lib -lm -ffast-math
	e-objcopy --srec-forceS3 --output-target srec e_raytrace.elf e_raytrace.srec

dump:
	e-objdump -d e_raytrace.elf

.PHONY: test
