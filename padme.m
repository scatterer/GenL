%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Gunnar K. Pálsson
% Description:      
%
% To do as user:    Nothing.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fpad = padme(F,thesize)

  half=thesize/2;
  Fpad                       = zeros(length(F)+half*2,1);
  Fpad(1:half)               = ones(half,1).*F(1);
  Fpad(half+1:end-half)      = F;
  Fpad(end-half+1:end)       = ones(half,1).*F(end);

end
