function fmm2dprini(unit1,unit2)
%FMM2DPRINI Initialize internal printing routines.
%
% Calling FMM2DPRINI(6,13) causes printing to screen and file fort.13.     
%

if (nargin == 1 )
unit2=0;
end

mex_id_ = 'prini(i int[x], i int[x])';
fmm2d(mex_id_, unit1, unit2, 1, 1);


