function [ier,box,center,corners]=d2tgetb(ibox,lists)
%D2TGETB Retrieve box information.
%
% [IER,BOX,CENTER,CORNERS]=D2TGETB(IBOX,LISTS);
%
% Input parameters:
%
% ibox - the box number for which the information is desired
% lists - storage area U.lists as created be D2TSTRCR or D2TSTRCREM
%
% Output parameters:
%
% ier - error return code
%    ier=0 - successful execution
%    ier=4 - ibox is either greater than the number of boxes 
%            in the structure or less than 1.
%
% box - an integer array dimensioned box(15). its elements describe 
%        the box number ibox, as follows:
%
%       1. level - the level of subdivision on which this box 
%             was constructed; 
%       2, 3  - the coordinates of this box among  all
%             boxes on this level
%       4 - the daddy of this box, identified by it address
%             in array boxes
%       5,6,7,8 - the  list of children of this box 
%             (eight of them, and the child is identified by its address
%             in the array boxes; if a box has only one child, only the
%             first of the four child entries is non-zero, etc.)
%       9 - the location in the array iz of the particles 
%             living in this box
%       10 - the number of particles living in this box
%       11 - the location in the array iztarg of the targets
%             living in this box
%       12 - the number of targets living in this box
%       13 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
%       14 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
%       15 - reserved for future use
%
% center - real (2) - the center of the box number ibox 
% corners - real (2,4) - the corners of the box number ibox 
%

ier = 0;
center = zeros(2,1);
corners = zeros(2,4);
box = zeros(1,15);

mex_id_ = 'd2tgetb(io int[x], i int[x], io int[], io double[], io double[], i double[])';
[ier, box, center, corners] = fmm2d(mex_id_, ier, ibox, box, center, corners, lists, 1, 1);



