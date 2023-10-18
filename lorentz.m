%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Anna L. Ravensburg, Gunnar K. Pálsson
% Description:      Calculate a Lorentzian to simulate the intensity originating from the substrate.
%
% To do as user:    Nothing.
%
% Note:
% Copyright (C) 2023 by the authors - All Rights Reserved
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details. A copy of the GNU
% General Public License can be obtained from the
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ip = lorentz(a,x)

	I  = a(1); % Intensity in cps
	w  = a(2); % peak width in degrees
	x0 = a(3); % interplanar spacing in substrate

    % calculate the final intensity
    Ip = I*1/pi*0.5*w./( (x-x0).^2 + (0.5*w).^2);

end
