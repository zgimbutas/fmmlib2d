%
%  Test the series expansion for Laplace and Helmholtz potentials in R^2
%
%  Currently, the Helmholtz FMM charge potential is positive near zero,
%  while the Laplace FMM charge potential is negative (we use log(r) for 
%  the Green's function). Same convention for the dipole potentials.
%
%
%> series(HankelH1(0,x),x);
%   -Pi + 2 I ln(2) - 2 I ln(x) - 2 I gamma
% - --------------------------------------- +
%                     Pi
%
%        -Pi + 2 I ln(2) - 2 I ln(x) + 2 I - 2 I gamma  2
%    1/4 --------------------------------------------- x  -
%                             Pi
%
%         -Pi + 2 I ln(2) - 2 I ln(x) + 3 I - 2 I gamma  4      6
%    1/64 --------------------------------------------- x  + O(x )
%                              Pi

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



euler_gamma=0.57721566490153286060651209008240243104

ima=1i;
helmholtz_charge_pot
laplace_charge_pot 
c=(-pi+2*ima*log(2)-2*ima*laplace_charge_pot-2*ima*euler_gamma)/pi
term0=-c *(ima/4);
error0=helmholtz_charge_pot-term0
term2=(c+2*ima/pi)/4 *(ima/4);
error2=helmholtz_charge_pot-term0-term2*h^2
term4=-(c+3*ima/pi)/64 *(ima/4);
error4=helmholtz_charge_pot-term0-term2*h^2-term4*h^4

