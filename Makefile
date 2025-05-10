OBJS = centroid_fold.o complex_number.o experiment.o fourier_transform.o main.o optimize.o mccaskill_1990_autodiff.o mccaskill_1990.o misc.o parameter.o real_logsumexp.o rintx.o sample_mccaskill.o test.o

CPPFLAGS = -I ./ -std=c++1y -fopenmp -static-libstdc++

.PHONY: release
release: EXE_FILE=rintc
release: CPPFLAGS+=-O2
release: build

.PHONY: debug
debug: EXE_FILE=rintc
debug: CPPFLAGS+=-g -O0
debug: build


.PHONY: build
build: ${OBJS}
	rm -f rintc debug
	g++ ${CPPFLAGS} ${OBJS} -o ${EXE_FILE}
	#rm -f *.o

.cpp.o:
	g++ ${CPPFLAGS} -c $< -o $@

clean:
	rm -f *.o rintc debug
