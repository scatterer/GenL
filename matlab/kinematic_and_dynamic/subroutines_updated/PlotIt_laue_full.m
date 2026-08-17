%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Anna L. Ravensburg, Gunnar K. Pálsson
% Description:      Plotting function used to plot fitting
%                   based on input "fit_laue_oscillations".
%
% To do as user:    Nothing.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotIt_laue_full(FVr_temp,iter,S_struct)

    global flag;

    %% SCATTERING PARAMETERS
    % FVr_temp: fitting parameters, S_struct: settings, variables provided
    
    Q = S_struct.instrument.Q; % in 1/Ångtröm
    
    theta    = S_struct.instrument.theta; % incident angle in Ångström
    twotheta = 2.*theta; % 2theta scattering angle

    %% CALCULATE MODEL
    [I,F_cost] = calc_model_full(FVr_temp,S_struct);
    disp_fitting_vector(FVr_temp,S_struct.stack,S_struct.instrument);
    %% PLOTTING
    if S_struct.GUI == 0 % access through command line version

        if S_struct.plot_unit == 0 % unit: 2Theta
            fig_fit = figure(1);
            set(fig_fit,'Position',[1,1,1200,800]);
            clf;
            set(fig_fit,'Color','w');
            plot(twotheta,S_struct.I(:),'Color','k'); % plot data
            hold on;
            plot(twotheta,I,'-r','LineWidth',1.5); % plot fit
        elseif S_struct.plot_unit == 1 % unit: Q
            plot(Q,S_struct.I(:),'Color','k'); % plot data
            hold on;
            plot(Q,I,'-r','LineWidth',1.5); % plot fit
        end
        set(gca,'YScale','log');
        xlabel('2Theta [degrees]'); % x axis lable
        ylabel('Intensity [cps]'); % y axis lable
        %title(strcat(S_struct.sample)); % title
        box on;
        set(gca,'FontSize',14);
        drawnow;
        pause(1e-9);

    elseif S_struct.GUI == 1 % access through GUI

        S_struct.axe.YLim = [-Inf Inf];
        if S_struct.plot_unit == 0 % unit: 2Theta
            cla(S_struct.axe);
            hold(S_struct.axe,'on');
            plot(S_struct.axe,twotheta,S_struct.I(:),'Color','k'); % plot data
            plot(S_struct.axe,twotheta,I, '-r','LineWidth',1.5); % plot fit
            hold(S_struct.axe,'off');
            S_struct.axe.XLabel.String = '2 Theta [degrees]';
        elseif S_struct.plot_unit == 1 % unit: Q
            cla(S_struct.axe);
            hold(S_struct.axe,'on');
            plot(S_struct.axe,Q,S_struct.I(:),'Color','k'); % plot data
            plot(S_struct.axe,Q,I, '-r','LineWidth',1.5); % plot fit
            hold(S_struct.axe,'off');
            S_struct.axe.XLabel.String = 'Q [1/Å]';
        end
        S_struct.axe.YScale = 'log';
        S_struct.axe.YLabel.String = 'Intensity [cps]';
    end
    
    
    %% SAVING
    if S_struct.saving == 1
        % save the plot
        %saveas(fig_fit, strcat(S_struct.path, S_struct.sample, '_LaueFit_graph'), 'png');
        
        % save fit
%        fileID = fopen(strcat(S_struct.path, S_struct.sample,'_LaueFit','.txt'),'w');
%        fprintf(fileID,strcat('Date: ', string(datetime("now")), '\n'));
%        fprintf(fileID,strcat('2Theta [degrees]     Q [1/Å]     I [cps]', '\n'));
%        for z = 1:length(twotheta)
%            fprintf(fileID,'%.7f    %.7f    %.7f', twotheta(z), Q(z), I(z));
%            fprintf(fileID,'\n');
%        end
%        fclose(fileID);
%    end

end