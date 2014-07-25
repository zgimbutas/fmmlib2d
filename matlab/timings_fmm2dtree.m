%
%  Test Helmholtz and Laplace particle FMMs in R^2
%

nsource = 100000

source = zeros(2,nsource);

  phi=rand(1,nsource)*2*pi;
  source(1,:)=.5*cos(phi);
  source(2,:)=.5*sin(phi);



%plot2(source(1,:),source(2,:))

%
%  build FMM tree structures in R^2
%
'build FMM tree in R^2'
tic
nbox = 100;
U = d2tstrcr(nsource,source,nbox);
time_fmmtree=toc


lused=U.lused
nboxes=U.nboxes
nlev=U.nlev
laddr=U.laddr(1:2,1:U.nlev)

