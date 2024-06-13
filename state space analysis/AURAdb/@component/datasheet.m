function UIFigure = datasheet(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%             ss = get(0,'screensize');     
          
            % Create UIFigure and hide until all components are created
            UIFigure = uifigure('Visible', 'off', 'Scrollable', 'on');
            UIFigure.Position = [100 100 640 751];
            UIFigure.Name = [obj.partNumber ' Datasheet'];
            
            %Find height for extra plots, if needed.
%             plotTypes = unique({obj.PlotData.Type});
            ExtraHeight = 300*floor((length(obj.graphs)-1)/2);

            % Create NameplateParams
            NameplateParams = uitable(UIFigure);
            NameplateParams.ColumnName = {'Parameter'; 'Value'; 'Unit'};
            NameplateParams.RowName = {};
            NameplateParams.Position = [40 445+ExtraHeight 547 121];
            
            
%             for i = 1:length(obj.parameters)
%                 
%             end
%             NameplateParams.Data = {'Rds(on) - Drain-Source On Resistance', obj.Ratings.ron*1000, ['m', char(937)]};


            % Create NameplateParametersLabel
            NameplateParametersLabel = uilabel(UIFigure);
            NameplateParametersLabel.Position = [78 565+ExtraHeight 129 22];
            NameplateParametersLabel.Text = 'Nameplate Parameters';

            P = properties(obj);
            P(strcmp(P, {'parameters'})) = [];
            P(strcmp(P, {'partNumber'})) = [];
            P(strcmp(P, {'graphs'})) = [];

            for i = 1:length(P)
                NameplateParams.Data{i,1} = P{i};
                NameplateParams.Data{i,2} = obj.(P{i});
            end
            


            % Create Characteristics
            Characteristics = uitable(UIFigure);
            Characteristics.ColumnName = {'Parameter'; 'Test Conditions'; 'Min'; 'Typ'; 'Max'; 'Unit'};
            Characteristics.RowName = {};
            Characteristics.Position = [40 287+ExtraHeight 547 121];
            
            Characteristics.Data = cell(length(obj.parameters),6);   
            
            for i = 1:length(obj.parameters)
                Characteristics.Data{i,1} = obj.parameters(i).name;
                Characteristics.Data{i,2} = char(obj.parameters(i).conditions(:));
                Characteristics.Data{i,3} = num2str(obj.parameters(i).min);
                Characteristics.Data{i,4} = num2str(obj.parameters(i).typ);
                Characteristics.Data{i,5} = num2str(obj.parameters(i).max);
                Characteristics.Data{i,6} = [obj.parameters(i).unit{1}, obj.parameters(i).unit{2}];
            end
%             Characteristics.Data = {'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
%                 'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
%                 'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];
%                 'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
%                 'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
%                 'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];
%                 'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
%                 'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
%                 'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];
%                 'Id', 'Continuous, T_A = 25C, Rtja = 220 C/W', '', 3.4, '', 'A'; ...
%                 'Id', 'Pulsed, T_A = 25C, Tpulse = 300 \mus', '', 28, '', 'A'; ... 
%                 'Rds(on)', 'Vgs = 5 V, Id = 1.5 A',24, 30, '', ['m', char(937)];};

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
     
            for i = 1:length(obj.graphs)
                Plot1 = uiaxes(UIFigure); 
%                 ind = find(strcmp({obj.PlotData.Type}, plotTypes{i}));
%                 title(Plot1, obj.PlotData(ind(1)).Title)
%                 xlabel(Plot1, obj.PlotData(ind(1)).Xlabel)
%                 ylabel(Plot1, obj.PlotData(ind(1)).Ylabel)
                obj.graphs(i).plot(Plot1);
                Plot1.Box = 'on';
                Plot1.Position = [40 + 300*mod(i-1,2), 50+ExtraHeight-250*floor((i-1)/2), 250, 200];
                
%                 for j = 1:length(ind)
%                     plot(Plot1, obj.PlotData(ind(j)).Data(:,1), obj.PlotData(ind(j)).Data(:,2));
%                     hold(Plot1, 'on');
% %                     legend(Plot1, {legend(),obj.PlotData(ind(j)).Content});
%                 end
%                 legend(Plot1, obj.PlotData(ind).Content);
            end

            % Create DatasheetLabel
            DatasheetLabel = uilabel(UIFigure);
            DatasheetLabel.Position = [40 710+ExtraHeight 547 22];
            DatasheetLabel.Text = obj.partNumber;
            DatasheetLabel.FontSize = 18;
            DatasheetLabel.FontWeight = 'bold';
            
            % Show the figure after all components are created
            scroll(UIFigure, DatasheetLabel.Position(1:2));
            UIFigure.Visible = 'on';
end

