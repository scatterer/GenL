%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Johan Bylin
% Description:      Get the coefficients for the form factor calculation.
%
% To do as user:    Nothing.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FF_coeff] = Read_form_factor_coefficients(Z, wavelength)

    function [f_1, f_2] = find_f1f2(Z, energy_source)
        
        % read in database BrennanCowanLong
        fid          = fopen( 'Form_Factor_and_Elemental_data/f1f2_BrennanCowanLong.dat');
        file_2       = textscan(fid,'%s','delimiter','\n');
        lines        = file_2{1};
        energy_dummy = 0;
        f_1_dummy    = 0;
        f_2_dummy    = 0;
        fclose(fid);
        
        % put the data in the correct format
        for j = 40:length(lines)
            if size(lines{j}) < 3
            else
                dummy_line = split(lines{j});
                if strcmp(dummy_line{1}, '#S')
                    if str2num(dummy_line{2}) == Z
                        
                        k = j + 5;
                        while size(lines{k}) > 0
                            dummy_line = split(lines{k});
                            energy_dummy  = [energy_dummy, str2double(dummy_line{1})];
                            f_1_dummy     = [f_1_dummy,    str2double(dummy_line{2})];
                            f_2_dummy     = [f_2_dummy,    str2double(dummy_line{3})];
                            
                            k = k + 1;
                        end
                        
                        [~, min_i] = min( abs(energy_dummy - energy_source));
                        f_1 = f_1_dummy(min_i);
                        f_2 = f_2_dummy(min_i);
                        break
                    end
                end
                
            end
        end
    end

% define necessary coefficients
hc            = 1.23984193*1e4; % h*c [eV*Å]
energy_source = hc/wavelength;

% create arrays
a   = zeros(5, length(Z));
b   = zeros(5, length(Z));
c   = zeros(1, length(Z));
f_1 = zeros(1, length(Z));
f_2 = zeros(1, length(Z));

% read in database Form_Factor_and_Elemental_data/ASF
file = importdata( 'Form_Factor_and_Elemental_data/ASF.DAT');

% find the form factor coefficients for the elements you are looking for
for i = 1:length(Z)
    
    a(:,i) = file.data(Z(i),1:5);
    b(:,i) = file.data(Z(i),7:11);
    c(:,i) = file.data(Z(i),6);
    
    [f_1(i), f_2(i)] = find_f1f2(Z(i), energy_source);
end

% return them
FF_coeff = struct('a', transpose(a), 'b', transpose(b), ...
    'c', transpose(c), 'f_1', transpose(f_1), ...
    'f_2', transpose(f_2));

end