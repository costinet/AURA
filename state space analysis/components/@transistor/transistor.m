classdef transistor < handle
    %transistor Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% Ratings
        partno
        Ratings
        
        %% Tables
        TableData
        
        %% Plots
        PlotData

    end
    
    methods
        function obj = transistor(partno, ron, qg, Rjc, Id, CossVds, IsdVsd, RonT, Vbr)          
            obj.partno = partno;
            obj.Ratings.ron = ron;
            obj.Ratings.qg = qg;
            obj.Ratings.Rjc = Rjc;
            obj.Ratings.Id = Id;
%             obj.PlotData(1) = CossVds;
%             obj.PlotData(2) = IsdVsd;
%             obj.PlotData(3) = RonT;% [RonT(:,1) RonT(:,2)*ron];
            obj.Ratings.Vbr = Vbr;
            
            Vds = linspace(0,Vbr,1000);
            
            

%             Coss = interp1(CossVds(:,1), CossVds(:,2), Vds);
%             CeqE = 2./Vds.^2.*cumtrapz(Vds, Coss.*Vds);
%             CeqE(1) = Coss(1);
%             CeqQ = 1./Vds.*cumtrapz(Vds, Coss);
%             CeqQ(1) = Coss(1);
%             obj.Vds = Vds;
%             obj.CeqQ = CeqQ;
%             obj.CeqE = CeqE;
        end
        
        function updateTemp(obj, T)
            obj.ron = min(interp1(obj.RonT(:,1), obj.RonT(:,2), T, 'linear', 'extrap'),obj.RonT(end,2));
        end
        
        function datasheet(obj)
%             ss = get(0,'screensize');     
%             [filepath,~,~] = fileparts(mfilename('fullpath'));
          
                 % Create UIFigure and hide until all components are created
            UIFigure = uifigure('Visible', 'off', 'Scrollable', 'on');
            UIFigure.Position = [100 100 640 751];
            UIFigure.Name = [obj.partno ' Datasheet'];
            
            %Find height for extra plots, if needed.
            plotTypes = unique({obj.PlotData.Type});
            ExtraHeight = 300*floor((length(plotTypes)-1)/2);

            % Create NameplateParams
            NameplateParams = uitable(UIFigure);
            NameplateParams.ColumnName = {'Parameter'; 'Value'; 'Unit'};
            NameplateParams.RowName = {};
            NameplateParams.Position = [40 445+ExtraHeight 547 121];
            
            NameplateParams.Data = {'Rds(on) - Drain-Source On Resistance', obj.Ratings.ron*1000, ['m', char(937)]};


            % Create NameplateParametersLabel
            NameplateParametersLabel = uilabel(UIFigure);
            NameplateParametersLabel.Position = [78 565+ExtraHeight 129 22];
            NameplateParametersLabel.Text = 'Nameplate Parameters';

            % Create Characteristics
            Characteristics = uitable(UIFigure);
            Characteristics.ColumnName = {'Parameter'; 'Test Conditions'; 'Min'; 'Typ'; 'Max'; 'Unit'};
            Characteristics.RowName = {};
            Characteristics.Position = [40 287+ExtraHeight 547 121];
            
            Characteristics.Data = {'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
                'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
                'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];
                'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
                'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
                'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];
                'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
                'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
                'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];
                'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
                'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
                'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];};

            % Create CharacteristicsLabel
            CharacteristicsLabel = uilabel(UIFigure);
            CharacteristicsLabel.Position = [78 407+ExtraHeight 85 22];
            CharacteristicsLabel.Text = 'Characteristics';
            
%             F = fieldnames(obj.PlotData);
%             types{1} = obj.PlotData.(F{1}).Type;
%             for i = 2:length(F)
%                 if isempty(obj.PlotData.(F{i}))
%                     types{i} = '';
%                     continue
%                 end
%                 types{i} = obj.PlotData.(F{i}).Type;
%             end
     
            for i = 1:length(plotTypes)
                Plot1 = uiaxes(UIFigure); 
                ind = find(strcmp({obj.PlotData.Type}, plotTypes{i}));
                title(Plot1, obj.PlotData(ind(1)).Title)
                xlabel(Plot1, obj.PlotData(ind(1)).Xlabel)
                ylabel(Plot1, obj.PlotData(ind(1)).Ylabel)   
                Plot1.Box = 'on';
                Plot1.Position = [40 + 300*mod(i-1,2), 50+ExtraHeight-250*floor((i-1)/2), 250, 200];
                
                for j = 1:length(ind)
                    plot(Plot1, obj.PlotData(ind(j)).Data(:,1), obj.PlotData(ind(j)).Data(:,2));
                    hold(Plot1, 'on');
%                     legend(Plot1, {legend(),obj.PlotData(ind(j)).Content});
                end
                legend(Plot1, obj.PlotData(ind).Content);
            end

            % Create DatasheetLabel
            DatasheetLabel = uilabel(UIFigure);
            DatasheetLabel.Position = [40 710+ExtraHeight 547 22];
            DatasheetLabel.Text = obj.partno;
            DatasheetLabel.FontSize = 18;
            DatasheetLabel.FontWeight = 'bold';
            
            % Show the figure after all components are created
            scroll(UIFigure, DatasheetLabel.Position(1:2));
            UIFigure.Visible = 'on';
        end
    end
    
end

