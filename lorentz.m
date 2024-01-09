%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Anna L. Ravensburg, Gunnar K. Pálsson
% Description:      Calculate a Lorentzian to simulate the intensity originating from the substrate.
%
% To do as user:    Nothing.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ip = lorentz(a,x)

	I  = a(1); % Intensity in cps
	w  = a(2); % peak width in degrees
	x0 = a(3); % interplanar spacing in substrate

    % calculate the final intensity
    Ip = I*1/pi*0.5*w./( (x-x0).^2 + (0.5*w).^2);

end
