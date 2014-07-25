%
%  Test Laplace particle target FMMs in R^2
%

nsource = 2000
ntarget = nsource*2

source = zeros(2,nsource);

  theta=rand(1,nsource)*pi;
  phi=rand(1,nsource)*2*pi;
  source(1,:)=.5*cos(phi);
  source(2,:)=.5*sin(phi);

target = zeros(2,ntarget);

  theta=rand(1,ntarget)*pi;
  phi=rand(1,ntarget)*2*pi;
  target(1,:)=.5*cos(phi) + 2;
  target(2,:)=.5*sin(phi);


%plot2(source(1,:),source(2,:))
%plot2(target(1,:),target(2,:))

%
%  Timings
%

ifcharge=1;
charge = ones(1,nsource)+1i*zeros(1,nsource);
ifdipole=1;
dipstr = ones(1,nsource)+1i*zeros(1,nsource);


ifcharge
ifdipole
ifpot = 1
ifgrad = 1
ifhess = 1
ifpottarg = 1
ifgradtarg = 1
ifhesstarg = 1


'Laplace particle target FMM in R^2'

tic
iprec=4
[U]=cfmm2dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
total_time=toc
speed=(nsource+ntarget)/total_time


'Laplace particle direct evaluation in R^2'

tic
[F]=c2dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
total_time=toc
speed=(nsource+ntarget)/total_time



if( ifpot ),
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
end

if( ifpot ),
re_rms_error_pot = norm(real(U.pot - F.pot),2)/sqrt(nsource)
end

if( ifgrad ),
rms_error_grad = norm((U.grad - F.grad),2)/sqrt(nsource)
end

if( ifhess ),
rms_error_hess = norm((U.hess - F.hess),2)/sqrt(nsource)
end



if( ifpottarg ),
rms_error_pottarg = norm((U.pottarg - F.pottarg),2)/sqrt(nsource)
end

if( ifpottarg ),
re_rms_error_pottarg = norm(real(U.pottarg - F.pottarg),2)/sqrt(nsource)
end

if( ifgradtarg ),
rms_error_gradtarg = norm((U.gradtarg - F.gradtarg),2)/sqrt(nsource)
end

if( ifhesstarg ),
rms_error_hesstarg = norm((U.hesstarg - F.hesstarg),2)/sqrt(nsource)
end



if( ifpot ),
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2)
end

if( ifpot ),
re_rel_error_pot = norm(real(U.pot - F.pot),2)/norm(real(F.pot),2)
end

if( ifgrad ),
rel_error_grad = norm((U.grad - F.grad),2)/norm((F.grad),2)
end

if( ifhess ),
rel_error_hess = norm((U.hess - F.hess),2)/norm((F.hess),2)
end



if( ifpottarg ),
rel_error_pottarg = norm((U.pottarg - F.pottarg),2)/norm((F.pottarg),2)
end

if( ifpottarg ),
re_rel_error_pottarg = norm(real(U.pottarg - F.pottarg),2)/norm(real(F.pottarg),2)
end

if( ifgradtarg ),
rel_error_gradtarg = norm((U.gradtarg - F.gradtarg),2)/norm((F.gradtarg),2)
end

if( ifhesstarg ),
rel_error_hesstarg = norm((U.hesstarg - F.hesstarg),2)/norm((F.hesstarg),2)
end


