# Helmholtz and Laplace FMM library in R^2.

Date: November 9, 2017

Version 1.2.2


The FMMLIB2D suite permits the evaluation of potential fields due to
particle sources, governed by either the Laplace or Helmholtz equation
in free space. The codes are easy to use and reasonably well optimized
for performance. A rudimentary manual is provided in the FMM2D/doc
directory.

FMMLIB2D contains both Fortran source code and versions compiled for
MATLAB under Mac OS X (64 bit), Windows (64 bit), and Linux (64 bit),
and Octave until Linux.


### License

```
Copyright (C) 2010-2012: Leslie Greengard and Zydrunas Gimbutas
Contact: greengard@cims.nyu.edu

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### Contents

```
src/ - Fortran source code
examples/ - Fortran testing drivers and makefiles
matlab/ - matlab scripts and mex files 
contrib/mwrap-0.33.3/ - mwrap source code
```

To test the library, please type `make test`, and when prompted
to `ENTER n` type something sensible like `10000` or `100000`.
You should see text output with small errors listed.
Warnings about floating-point exceptions are normal and to be ignored.

**Note** if using gfortran v 10 or above: Since we use passing of size-1 arrays
as pointers, GCC10+ raises errors. You will need to add
```
FFLAGS+=-std=legacy
```
in the relevant sections of `src/Makefile` and `examples/*.make`
which turns these into mere warnings.



### Fortran

Particle FMM routines.

```
hfmm2dpartself - Helmholtz particle FMM in R^2.
lfmm2dpartself - Laplace particle FMM in R^2 (complex).
rfmm2dpartself - Laplace particle FMM in R^2 (real).
zfmm2dpartself - Laplace particle FMM in R^2 (Cauchy sums).
cfmm2dpartself - Laplace particle FMM in R^2 (generalized Cauchy sums).

hfmm2dparttarg - Helmholtz particle target FMM in R^2.
lfmm2dparttarg - Laplace particle target FMM in R^2 (complex).
rfmm2dparttarg - Laplace particle target FMM in R^2 (real).
zfmm2dparttarg - Laplace particle target FMM in R^2 (Cauchy sums).
cfmm2dparttarg - Laplace particle target FMM in R^2 (generalized Cauchy sums).
```

Direct evaluation routines.

```
h2dpartdirect - Helmholtz particle target interactions in R^2.
l2dpartdirect - Laplace particle target interactions in R^2 (complex).
r2dpartdirect - Laplace particle target interactions in R^2 (real).
z2dpartdirect - Laplace particle target interactions in R^2 (Cauchy).
c2dpartdirect - Laplace particle target interactions in R^2 (g. Cauchy).
```

### Matlab

```
% Helmholtz and Laplace FMMs in R^2.
%
% Particle FMM routines.
%   hfmm2dpart      - Helmholtz particle FMM in R^2.
%   lfmm2dpart      - Laplace particle FMM in R^2 (complex).
%   rfmm2dpart      - Laplace particle FMM in R^2 (real).
%   zfmm2dpart      - Laplace particle FMM in R^2 (Cauchy sums).
%   cfmm2dpart      - Laplace particle FMM in R^2 (g. Cauchy sums).
%
% Direct evaluation routines.
%   h2dpartdirect - Helmholtz particle target interactions in R^2.
%   l2dpartdirect - Laplace particle target interactions in R^2 (complex).
%   r2dpartdirect - Laplace particle target interactions in R^2 (real).
%   z2dpartdirect - Laplace particle target interactions in R^2 (Cauchy).
%   c2dpartdirect - Laplace particle target interactions in R^2 (g. Cauchy).
%
```

### Acknowledgments

This work was supported in part by the Department of Energy and in
part by the Air Force Office of Scientific Research under MURI grant
FA9550-06-1-0337 and NSSEFF Program Award FA9550-10-0180, in part by
the NSF under grant DMS09-34733, and in part by Meyer Sound
Laboratories, Inc.

