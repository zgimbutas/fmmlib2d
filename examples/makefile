
hfmm2dpart:
	make -f hfmm2dpart.make -j4

cfmm2dpart:
	make -f cfmm2dpart.make -j4

zfmm2dpart:
	make -f zfmm2dpart.make -j4

lfmm2dpart:
	make -f lfmm2dpart.make -j4

rfmm2dpart:
	make -f rfmm2dpart.make -j4

hfmm2dpart-openmp:
	make -f hfmm2dpart.make HOST=linux-gfortran-openmp -j4

cfmm2dpart-openmp:
	make -f cfmm2dpart.make HOST=linux-gfortran-openmp -j4

zfmm2dpart-openmp:
	make -f zfmm2dpart.make HOST=linux-gfortran-openmp -j4

lfmm2dpart-openmp:
	make -f lfmm2dpart.make HOST=linux-gfortran-openmp -j4

rfmm2dpart-openmp:
	make -f rfmm2dpart.make HOST=linux-gfortran-openmp -j4

clean: 
	rm -f *.o
 
distclean: 
	rm -f int2 fort.* *.o

