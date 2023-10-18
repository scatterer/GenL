%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:         I_z = left_win(S_x,S_y)
% Author:           Rainer Storn
% Description:      left_win(S_x,S_y) takes structures S_x and S_y as an argument.
%                   The function returns 1 if the left structure of the input structures,
%                   i.e. S_x, wins. If the right structure, S_y, wins, the result is 0.
% Parameters:       S_x.I_nc     (I)    Number of constraints (must be the same for x and y).
%                   S_x.I_no     (I)    Number of objectives (must be the same for x and y).
%                   S_x.FVr_ca   (I)    Constraint array containing the constraint violation values.
%                                       If the value is 0 the constraint is met. If it is > 0 it is
%                                       still violated.
%                   S_x.FVr_oa   (I)    Objective array containing cost values which are supposed to be
%                                       minimized.
% Return value:     I_z          (O)    If S_x wins over S_y then I_z=1 else I_z=0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 1, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. A copy of the GNU 
% General Public License can be obtained from the 
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_z = left_win(S_x,S_y);
I_z = 1;  %start with I_z=1

%----deal with the constraints first. If constraints are not met------
%----S_x can't win.---------------------------------------------------
if (S_x.I_nc > 0)
   for k=1:S_x.I_nc
      if (S_x.FVr_ca(k) > 0) %if constraint is not yet met
        if (S_x.FVr_ca(k) > S_y.FVr_ca(k))%if just one constraint of S_x is not improved
           I_z = 0;
        end
      end
   end   
end

if (S_x.I_no > 0)
   for k=1:S_x.I_no
      if (S_x.FVr_oa(k) > S_y.FVr_oa(k))%if just one objective of S_x is less
         I_z = 0;
      end
   end
end


