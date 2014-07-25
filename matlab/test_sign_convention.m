%
%  Test the sign conventions for Laplace and Helmholtz FMM in R^2
%
%  Currently, the Helmholtz FMM charge potential is positive near zero,
%  while the Laplace FMM charge potential is negative (we use log(r) for 
%  the Green's function). Same convention for the dipole potentials.
%


nsource = 1
ntarget = 4

source = zeros(2,nsource);
target = zeros(2,ntarget);

phi=(0:ntarget-1)/ntarget*pi*2;

h=1e-3;
target(1,:)=h*cos(phi);
target(2,:)=h*sin(phi);

source
target


ifcharge=1;
charge = ones(1,nsource)+1i*zeros(1,nsource);
ifdipole=0;
dipstr = ones(1,nsource)+1i*zeros(1,nsource);
dipvec = ones(2,nsource);

ifpot = 0;
ifgrad = 0;
ifhess = 0;
ifpottarg = 1;
ifgradtarg = 0;
ifhesstarg = 0;

%  check charge convention

zk=1

tic
[F]=h2dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);

helmholtz_charge_pot=F.pottarg

tic
[F]=l2dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);

laplace_charge_pot=F.pottarg

%  check dipole convention

ifcharge=0;
charge = ones(1,nsource)+1i*zeros(1,nsource);
ifdipole=1;
dipstr = ones(1,nsource)+1i*zeros(1,nsource);
dipvec = ones(2,nsource);

zk=1

tic
[F]=h2dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);

helmholtz_dipole_pot=F.pottarg

tic
[F]=l2dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);

laplace_dipole_pot=F.pottarg

