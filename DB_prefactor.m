%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Anna L. Ravensburg, Gunnar K. Pálsson
% Description:      Get the element specific Debye-Waller pre-factor.
%                   
%
% To do as user:    Nothing.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dwpf = DB_prefactor(Z)

    % set dwpf to 0 as a default for it not working
    dwpf = 0;

    % if element is in data base, but it does not work for alloys
    if length(Z) == 1
        if Z == 3 % Li bcc
            dwpf = -4.9;
        elseif Z == 11 % Na bcc
            dwpf = -6.6;
        elseif Z == 13 % Al fcc
            dwpf = -0.81;
        elseif Z == 14 % Si diamond
            dwpf = -0.52;
        elseif Z == 19 % K bcc
            dwpf = -10.7;
        elseif Z == 23 % V bcc
            dwpf = -0.58;
        elseif Z == 24 % Cr bcc
            dwpf = -0.263;
        elseif Z == 26 % Fe bcc
            dwpf = -0.34;
        elseif Z == 28 % Ni fcc
            dwpf = -0.34;
        elseif Z == 29 % Cu fcc
            dwpf = -0.56;
        elseif Z == 32 % Ge diamond
            dwpf = -0.61;
        elseif Z == 41 % Nb bcc
            dwpf = -0.45;
        elseif Z == 42 % Mo bcc
            dwpf = -0.220;
        elseif Z == 46 % Pd fcc
            dwpf = -0.45;
        elseif Z == 47 % Ag fcc
            dwpf = -0.73;
        elseif Z == 73 % Ta bcc
            dwpf = -0.32;
        elseif Z == 74 % W bcc
            dwpf = -0.161;
        elseif Z == 78 % Pt fcc
            dwpf = -0.37;
        elseif Z == 79 % Au fcc
            dwpf = -0.62;
        elseif Z == 82 % Pb fcc
            dwpf = -2.11;    
        end       
    end
end