%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Anna L. Ravensburg, Gunnar K. Pálsson
% Description:      Read poscar for atomic positions of the substrate.
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
function [type,type_nr,r,a1,a2,a3] = read_poscar(filename)
    filename = strcat('POSCAR/',filename);
    structure = importdata(filename,' ',8);
    scaling   = str2num(structure.textdata{2});
    a1        = str2num(structure.textdata{3})*scaling;
    a2        = str2num(structure.textdata{4})*scaling;
    a3        = str2num(structure.textdata{5})*scaling;
    typedata    = structure.textdata;

    % Elements, second line from the bottom of header
    types   = typedata{end-2}; type = split(types); 

    % How many elements of each one line from the bottom of header
    type_nr = typedata{end-1};
    type_nr = str2num(type_nr); % Get the number of each species
    r       = structure.data;
    r = double(r);
end