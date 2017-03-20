function [U]=r2dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg)
%R2DPARTDIRECT Laplace interactions in R^2, direct evaluation (real).
%
% Laplace direct evaluation in R^2: evaluate all pairwise particle
% interactions (ignoring self-interactions) and interactions with targets.
%
% r2d: charge and dipstr are real valued, x \in R^2
%
% \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
%                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
%
% or, more precisely,
%
% \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
%                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
%
%
% [U]=R2DPARTDIRECT(NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC);
%
% [U]=R2DPARTDIRECT(NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFGRAD,IFHESS);
%
% [U]=R2DPARTDIRECT(NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFGRAD,IFHESS,...
%         NTARGET,TARGET,IFPOTTARG,IFGRADTARG,IFHESSTARG);
%
%
% This subroutine evaluates the Laplace potential and gradient due
% to a collection of charges and dipoles. We use log(r) for the 
% Green's function. Self-interactions are not-included.
%
% Input parameters:
% 
% nsource - number of sources
% source - real (2,nsource): source locations
% ifcharge - charge computation flag
%
%         0 => do not compute
%         1 => include charge contribution
% 
% charge - real (nsource): charge strengths 
% ifdipole - dipole computation flag
%
%         0 => do not compute
%         1 => include dipole contributions
% 
% dipole - real (nsource): dipole strengths
% dipvec - real (2,source): dipole orientation vectors
%
% ifpot - potential computation flag, 1 => compute the potential, otherwise no
% ifgrad - gradient computation flag, 1 => compute the gradient, otherwise no
% ifhess - hessian computation flag, 1 => compute the hessian, otherwise no
%
% ntarget - number of targets
% target - real (2,ntarget): target locations
%
% ifpottarg - target potential computation flag, 
%      1 => compute the target potential, otherwise no
% ifgradtarg - target gradient computation flag, 
%      1 => compute the target gradient, otherwise no
% ihesstarg - target hessian computation flag 
%      1 => compute the hessian, otherwise no
%
% Output parameters: 
%
% U.pot - real (nsource) - potential at source locations
% U.grad - real (2,nsource) - gradient  at source locations
% U.hess - real (3,nsource) - hessian at source locations
% U.pottarg - real (ntarget) - potential at target locations
% U.gradtarg - real (2,ntarget) - gradient  at target locations
% U.hesstarg - real (3,ntarget) - hessian at target locations
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%

if( nargin == 7 ) 
  ifpot = 1;
  ifgrad = 1;
  ifhess = 1;
  ntarget = 0;
  target = zeros(2,1);
  ifpottarg = 0;
  ifgradtarg = 0;
  ifhesstarg = 0;
end

if( nargin == 10 ) 
  ntarget = 0;
  target = zeros(2,1);
  ifpottarg = 0;
  ifgradtarg = 0;
  ifhesstarg = 0;
end

if( nargin == 12 ) 
  ifpottarg = 1;
  ifgradtarg = 1;
  ifhesstarg = 1;
end


ifcharge = double(ifcharge); ifdipole = double(ifdipole);
ifpot = double(ifpot); ifgrad = double(ifgrad); ifhess = double(ifhess);
ifpottarg = double(ifpottarg); ifgradtarg = double(ifgradtarg);
ifhesstarg = double(ifhesstarg);

pot=0;
grad=zeros(2,1);
hess=zeros(3,1);
pottarg=0;
gradtarg=zeros(2,1);
hesstarg=zeros(3,1);

if( ifpot == 1 ), pot=zeros(1,nsource); end;
if( ifgrad == 1 ), grad=zeros(2,nsource); end;
if( ifhess == 1 ), hess=zeros(3,nsource); end;
if( ifpottarg == 1 ), pottarg=zeros(1,ntarget); end;
if( ifgradtarg == 1 ), gradtarg=zeros(2,ntarget); end;
if( ifhesstarg == 1 ), hesstarg=zeros(3,ntarget); end;

ier=0;

mex_id_ = 'r2dpartdirect(i int[x], i double[xx], i int[x], i double[], i int[x], i double[], i double[xx], i int[x], io double[], i int[x], io double[], i int[x], io double[], i int[x], i double[], i int[x], io double[], i int[x], io double[], i int[x], io double[])';
[pot, grad, hess, pottarg, gradtarg, hesstarg] = fmm2d(mex_id_, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, ifgrad, grad, ifhess, hess, ntarget, target, ifpottarg, pottarg, ifgradtarg, gradtarg, ifhesstarg, hesstarg, 1, 2, nsource, 1, 1, 2, nsource, 1, 1, 1, 1, 1, 1, 1);


if( ifpot == 1 ) U.pot=pot; end
if( ifgrad == 1 ) U.grad=grad; end
if( ifhess == 1 ) U.hess=hess; end
if( ifpottarg == 1 ) U.pottarg=pottarg; end
if( ifgradtarg == 1 ) U.gradtarg=gradtarg; end
if( ifhesstarg == 1 ) U.hesstarg=hesstarg; end
U.ier=ier;



