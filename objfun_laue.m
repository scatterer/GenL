%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Anna L. Ravensburg, Gunnar K. Pálsson
% Description:      Objective function used to find optimal fitting
%                   parameters based on input "fit_laue_oscillations".
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

function S_MSE= objfun_laue(FVr_temp, S_struct)

    %% SCATTERING PARAMETERS
    % FVr_temp are the fitting parameters, S_struct are the parameter provided
    % to this program.
    
    Q = S_struct.Q; % in 1/Ångtröm
    f = S_struct.f;
    
    % debye-waller factor
    DebyeWaller = 1*exp(S_struct.dwf*(Q/4/pi).^2);
    % intensity
    I = zeros(size(Q));
    
    %% SAMPLE PARAMETERS
    d     = FVr_temp(1); % interplanar spacing in Ångström

    % strain
    alpha_1 = FVr_temp(2); % 1st strain parameter
    alpha_2 = FVr_temp(8); % 2nd strain parameter
    
    % background
    a     = FVr_temp(3); % 1st background parameter
    b     = FVr_temp(4); % 2nd background parameter
    if S_struct.background == 1
        c     = FVr_temp(9); % 2nd background parameter for polynomial bkg
    end
    
    % layer roughness
    sigma = FVr_temp(10); % in Ångström

    % rounded number of layers
    Nlayers = round(FVr_temp(5));

    % calculate DWF
    DF = DebyeWaller.*f;
    
    %% STRUCTURE FACTOR
    % if layer roughness should be included in the fitting
    if S_struct.roughness == 1 && sigma > 0

        % range of 3 sigma plus/minus for roughness to be included
        nvector = max(1,round(Nlayers-3*sigma)):1:max(1,round(Nlayers+3*sigma));

        % out-of-plane positions of the atoms
        x   = cell(size(nvector));
        % positions
        pos = cell(size(nvector));
        % including roughness
        B   = zeros(size(nvector));
        A   = 1/(sqrt(2*pi)*sigma);

        % strain
        strain   = cell(size(nvector));


        for N = nvector
        
            % out-of-plane positions of the atoms
            x{N} = ones(1,N)*d;
            % position
            p      = cumsum(x{N});
            pos{N} = [0 p];
            % start position
            pos{N} = [0 p];

            % including roughness
            B(N) = A*exp(-(N-Nlayers)^2/(2*sigma^2));
        
            % strain
            if S_struct.strained == 1
                % change by strain
                strain{N} = exp(-alpha_1*pos{N})+exp(-alpha_2*pos{N});
                % changes positions by strain
                pos{N} = pos{N} - strain{N};
            elseif S_struct.strained == -1
                % change by strain
                strain{N} = exp(-alpha_1*pos{N})+exp(-alpha_2*pos{N});
                % changes positions by strain
                pos{N} = pos{N} + strain{N};
            end
        end
        
        % initialization structure factor
        F = zeros(size(Q));
    
        % calculate structure factor
        for N = nvector
            for i =1:length(Q)
                F(i) = sum(exp(1i*Q(i).*pos{N}));
            end
            I = I + abs(F).^2*B(N); % calculated Intensity
        end

        % calculate intensity
        I = (DF).^2.*I;

    % if layer roughness should not be included in the fitting    
    elseif S_struct.roughness == 0

        % out-of-plane positions of the atoms
        x = ones(1,round(FVr_temp(5)))*d;
        pos = cumsum(x);
        % start position
        pos = [0 pos];
        
        % strain
        if S_struct.strained == 1
            % change by strain
            strain = exp(-alpha_1*pos)+exp(-alpha_2*pos);
            % changes positions by strain
            pos = pos - strain;
        elseif S_struct.strained == -1
            % change by strain
            strain = exp(-alpha_1*pos)+exp(-alpha_2*pos);
            % changes positions by strain
            pos = pos + strain;
        end

        % initialization structure factor
        F = zeros(size(Q));

        % calculate structure factor
        for i =1:length(Q)
            F(i) = DF(i)*sum(exp(1i*Q(i).*pos));
        end

        I = abs(F).^2; % calculated Intensity
    end
    
    %% STACK THE ATOMS FINAL
    % out-of-plane positions of the atoms
    x = ones(1,round(FVr_temp(5)))*d;
    pos = cumsum(x);
    % start position
    pos = [0 pos];
    
    % strain
    if S_struct.strained == 1
        % change by strain
        strain = exp(-alpha_1*pos)+exp(-alpha_2*pos);
        % changes positions by strain
        pos = pos - strain;
    elseif S_struct.strained == -1
        % change by strain
        strain = exp(-alpha_1*pos)+exp(-alpha_2*pos);
        % changes positions by strain
        pos = pos + strain;
    end
    
    %% CALCULATE SCATTERING PROPERTIES
    theta = asind(Q/4/pi*S_struct.lambda); % incident angle in Ångström
    twotheta = 2.*theta;
    theta_m = S_struct.theta_m;
    tau = pos(end)*1e-10; % film thickness in m
    mu = S_struct.mu; % absorption coefficient
    
    %% CALCULATE THE INTENSITY

    % optics on incident side
    if S_struct.theta_mPath == 1
        % with linear background a*Q+b
        if S_struct.background == 0
            I    = FVr_temp(7)*GaussConv(Q,I,FVr_temp(6)).* (1-exp(-2*mu*tau./sind(theta))).* (1+cosd(2*theta_m)^2.*cosd(2*theta).^2)./(1+cosd(2*theta_m).^2).* 1./(sind(2*theta)) + a*Q+b;
        % with polynomial background c*(a+Q).^2+b
        elseif S_struct.background == 1
            I    = FVr_temp(7)*GaussConv(Q,I,FVr_temp(6)).* (1-exp(-2*mu*tau./sind(theta))).* (1+cosd(2*theta_m)^2.*cosd(2*theta).^2)./(1+cosd(2*theta_m).^2).* 1./(sind(2*theta)) + c*(a+Q).^2+b;
        end

    % optics on detector side 
    elseif S_struct.theta_mPath == -1
        % with linear background a*Q+b
        if S_struct.background == 0
            I    = FVr_temp(7)*GaussConv(Q,I,FVr_temp(6)).* (1-exp(-2*mu*tau./sind(theta))).* 0.5*(1+cosd(2*theta).^2*cosd(2*theta_m).^2).* 1./(sind(2*theta)) + a*Q+b;
        % with polynomial background c*(a+Q).^2+b
        elseif S_struct.background == 1
            I    = FVr_temp(7)*GaussConv(Q,I,FVr_temp(6)).* (1-exp(-2*mu*tau./sind(theta))).* 0.5*(1+cosd(2*theta).^2*cosd(2*theta_m).^2).* 1./(sind(2*theta)) + c*(a+Q).^2+b;
        end
    end
    
    %% INCLUDING SCATTERING FROM THE SUBSTRATE
    if S_struct.substrate == 1
        a_lor     = [FVr_temp(11) FVr_temp(12) FVr_temp(13)];
        I_lor     = lorentz(a_lor,Q);
        I         = I_lor + I;
    end
        
    %% FIT FUNCTION
    F_cost = sum( abs( log10(I(:)) - log10(S_struct.I) ))./(S_struct.N-1);
    
    % cost function settings
    S_MSE.I_nc      = 0; % no constraints
    S_MSE.FVr_ca    = 0; % no constraint array
    S_MSE.I_no      = 1; % number of objectives (costs)
    S_MSE.FVr_oa(1) = F_cost;
end