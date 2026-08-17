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
% Sears Acta Cryst. A47, 441 (1991)
% Room temperature DW factors of elemental crystals

% DW = exp(-2*W) (is the intensity)
% 2*W = (q*u)^2 , q = 4*pi/lambda*sind(theta)
% 2*W = 2*B*[sind(theta)/lambda]^2
% B   = 8*pi^2*u^2

% W in Warren's book is called M_n where n labels the atom in the unit cell
% We replace the form factor of each element f_n by f_n*exp(-M_n) according to Warren

function B = DB_prefactor(Z)

    % set dwpf to 0 as a default for it not working
    B = 0;

    % if element is in data base, but it does not work for alloys
    % Calculated values by Sears:
    % if length(Z) == 1
    %     if Z == 3 % Li bcc
    %         dwpf = -4.9;
    %     elseif Z == 11 % Na bcc
    %         dwpf = -6.6;
    %     elseif Z == 13 % Al fcc
    %         dwpf = -0.81;
    %     elseif Z == 14 % Si diamond
    %         dwpf = -0.52;
    %     elseif Z == 19 % K bcc
    %         dwpf = -10.7;
    %     elseif Z == 23 % V bcc
    %         dwpf = -0.58;
    %     elseif Z == 24 % Cr bcc
    %         dwpf = -0.263;
    %     elseif Z == 26 % Fe bcc
    %         dwpf = -0.34;
    %     elseif Z == 28 % Ni fcc
    %         dwpf = -0.34;
    %     elseif Z == 29 % Cu fcc
    %         dwpf = -0.56;
    %     elseif Z == 32 % Ge diamond
    %         dwpf = -0.61;
    %     elseif Z == 41 % Nb bcc
    %         dwpf = -0.45;
    %     elseif Z == 42 % Mo bcc
    %         dwpf = -0.220;
    %     elseif Z == 46 % Pd fcc
    %         dwpf = -0.45;
    %     elseif Z == 47 % Ag fcc
    %         dwpf = -0.73;
    %     elseif Z == 73 % Ta bcc
    %         dwpf = -0.32;
    %     elseif Z == 74 % W bcc
    %         dwpf = -0.161;
    %     elseif Z == 78 % Pt fcc
    %         dwpf = -0.37;
    %     elseif Z == 79 % Au fcc
    %         dwpf = -0.62;
    %     elseif Z == 82 % Pb fcc
    %         dwpf = -2.11;    
    %     end       
    % end
    % Experimental values from Sears paper
    if length(Z) == 1
      if Z == 3 % Li bcc
        B  = 4.1;
      elseif Z == 11 % Na bcc
        B  = 7.9;
      elseif Z == 13 % Al fcc
        B  = 0.86;
      elseif Z == 14 % Si diamond
        B = 0.45;
      elseif Z == 19 % K bcc
        B = 12;
      elseif Z == 23 % V bcc
        B = 0.55;
      elseif Z == 24 % Cr bcc
        B = 0.26;
      elseif Z == 26 % Fe bcc
        B = 0.35;
      elseif Z == 28 % Ni fcc
        B = 0.37;
      elseif Z == 29 % Cu fcc
        B = 0.57;
      elseif Z == 32 % Ge diamond
        B = 0.57;
      elseif Z == 41 % Nb bcc
        B = 0.49;
      elseif Z == 42 % Mo bcc
        B = 0.25;
      elseif Z == 46 % Pd fcc
        B = 0.45;
      elseif Z == 47 % Ag fcc
        B = 0.79;
      elseif Z == 73 % Ta bcc
        B = 0.32;
      elseif Z == 74 % W bcc
        B = 0.18;
      elseif Z == 78 % Pt fcc
        B = 0.32;
      elseif Z == 79 % Au fcc
        B = 0.57;
      elseif Z == 82 % Pb fcc
        B = 2.42;
      end
    end
end