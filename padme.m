%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Gunnar K. Pálsson
% Description:      
%
% To do as user:    Nothing.
%
% Note:
% Copyright (C) 2023 by the author - All Rights Reserved
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details. A copy of the GNU
% General Public License can be obtained from the
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fpad = padme(F,thesize)

  half=thesize/2;
  Fpad                       = zeros(length(F)+half*2,1);
  Fpad(1:half)               = ones(half,1).*F(1);
  Fpad(half+1:end-half)      = F;
  Fpad(end-half+1:end)       = ones(half,1).*F(end);

end
