%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Gunnar K. Pálsson
% Description:      Calculate a Gaussian convolution for the final
%                   scattering intensity.
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

function Iconved = GaussConv(X,Y,fwhm)
  
  sigma  = fwhm/2/sqrt(2*log(2));
  points = floor(4*sigma./abs(X(1)-X(2)));

  thesize = 154;
  
  Fpad    = padme(Y,thesize);
    
  dtth   = (-points:points+1)*(abs(X(1)-X(2)));

  G      = 1./(sigma*sqrt(2*pi)).*exp(-(dtth./(sqrt(2)*sigma)).^2);
  G      = G./sum(G);

  Ipad    = conv(Fpad,G,'same');  
  Iconved = Ipad(thesize/2+1:end-thesize/2);
end
