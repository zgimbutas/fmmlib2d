%
%  Test Helmholtz particle target FMMs in R^2
%

nsource = 100000
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

zk=1;
ifcharge=1;
charge = ones(1,nsource)+1i*zeros(1,nsource);
ifdipole=1;
dipstr = ones(1,nsource)+1i*zeros(1,nsource);
dipvec = ones(2,nsource);


zk
ifcharge
ifdipole
ifpot = 1
ifgrad = 1
ifhess = 1
ifpottarg = 1
ifgradtarg = 1
ifhesstarg = 1

'Helmholtz particle FMM in R^2'

tic
iprec=4
[U]=hfmm2dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess);
total_time=toc
speed=(nsource)/total_time

'Helmholtz particle target FMM in R^2'

tic
iprec=4
[U]=hfmm2dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
total_time=toc
speed=(nsource+ntarget)/total_time

