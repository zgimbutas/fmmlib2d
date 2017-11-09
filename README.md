This is the second public release of the particle FMM library in R^2.

Date: November 8, 2017

Version 1.2.1

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

3. Neither the name of the copyright holder nor the names of its
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

To test the library, please type `make test`. 


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
