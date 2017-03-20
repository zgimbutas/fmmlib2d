function [ier,list,nlist]=d2tgetl(ibox,itype,lists)
%D2TGETL Retrieve list information.
%
% [IER,LIST,NLIST]=D2TGETL(IBOX,ITYPE,LISTS);
%
% Input parameters:
%
% ibox - the box number for which the information is desired
% itype - the type of the desired list for the box ibox
% lists - storage area U.lists as created be D2TSTRCR or D2TSTRCREM
%
% Output parameters:
%
% ier - the error return code.
%    ier=0 - successful execution
%    ier=4 - the list  itype  for the box  ibox  is empty
% list - the list  itype  for the box  ibox 
% nlist - the number of elements in array  list
%

ier = 0;
list = zeros(1,10000);
nlist = 0;

mex_id_ = 'd2tgetl(io int[x], i int[x], i int[x], io int[], io int[], i double[])';
[ier, list, nlist] = fmm2d(mex_id_, ier, ibox, itype, list, nlist, lists, 1, 1, 1);

if( ier == 0 ) list = list(1,1:nlist); end
if( ier  > 0 ) list = list(1,1); end

