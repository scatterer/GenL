%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Johan Bylin
% Description:      Get the Q dependent form factors for a certain material out, if you
%                   know the elements (and their Z) and their composition in the material.
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
function [f, f_sqrd_real, f_av_sqrd_real] = Form_factors(Q, FF, composition)

% Q dependent form factors
form_factor = zeros(length(composition), length(Q));
f           = zeros(length(Q), 1);
f_sqrd      = zeros(length(Q), 1);

% calculate contributions based on composition
for i = 1:length(composition)
    form_factor(i,:) = FF.c(i) + FF.f_1(i) + 1j*FF.f_2(i); 
    for j = 1:5
        gaussian_basis = FF.a(i,j).*exp( -FF.b(i,j).*(Q/(4.*pi)).^2);
        form_factor(i,:) = form_factor(i,:) + transpose(gaussian_basis);
    end
end

for i=1:length(composition)
    f      = f + transpose( composition(i).*form_factor(i,:));
    f_sqrd = f_sqrd + transpose(composition(i).*form_factor(i,:).*conj(form_factor(i,:)) );
end

% get the real parts of the form factors
f_sqrd_real    = real(f_sqrd);
f_av_sqrd_real = real(f.*conj(f));


