function UIFigure = datasheet(obj)
%datasheet Creates a datasheet of the component

%             ss = get(0,'screensize');
try
    % Create UIFigure and hide until all components are created
    UIFigure = uifigure('Visible', 'off', 'Scrollable', 'on');
    UIFigure.Position = [100 100 640 751];
    UIFigure.Name = [obj.partNumber ' Datasheet'];

    %Find height for extra plots, if needed.
    %             plotTypes = unique({obj.PlotData.Type});
    ExtraHeight = 300*floor((length(obj.graphs)-1)/2);


    % Create AURA Label
    AURA = uilabel(UIFigure);
    AURA.HorizontalAlignment = 'center';
    AURA.FontName = 'Courier';
    AURA.FontSize = 18;
    AURA.Position = [221 590+ExtraHeight 337 132];
    AURA.Text = {'     _   _   _  ____    _    '; '    / \ | | | |/ _  |  / \   '; '   / _ \| | | | (_| | / _ \  '; '  / ___ | |_| |> _  |/ ___ \ '; ' /_/   \_\___//_/ |_/_/   \_\'};

    warning ('off','MATLAB:ui:Image:invalidIconNotInPath');
    Image = uiimage(UIFigure);
    Image.Position = [71 590+ExtraHeight 100 100];
    Image.ImageSource = [obj.partNumber '.png'];

    % Create NameplateParams
    NameplateParams = uitable(UIFigure);
    NameplateParams.ColumnName = {'Parameter'; 'Value'; 'Unit'};
    NameplateParams.RowName = {};
    NameplateParams.Position = [40 445+ExtraHeight 547 121];

    %            if isequal(class(obj),'transistor')
    %                NameplateParams.Data = {'Rds(on) - Drain-Source On Resistance', obj.Ratings.ron*1000, ['m', char(937)]};
    %            end
    %             for i = 1:length(obj.parameters)
    %
    %             end
    %             NameplateParams.Data = {'Rds(on) - Drain-Source On Resistance', obj.Ratings.ron*1000, ['m', char(937)]};

    %if isequal(class(obj),'transistor')




    % Create NameplateParametersLabel
    NameplateParametersLabel = uilabel(UIFigure);
    NameplateParametersLabel.Position = [78 565+ExtraHeight 129 22];
    NameplateParametersLabel.Text = 'Nameplate Parameters';

    % Create Characteristics
    Characteristics = uitable(UIFigure);
    Characteristics.ColumnName = {'Parameter'; 'Test Conditions'; 'Min'; 'Typ'; 'Max'; 'Unit'};
    Characteristics.RowName = {};
    Characteristics.Position = [40 287+ExtraHeight 547 121];

    Characteristics.Data = cell(length(obj.parameters),6);
    if isequal(class(obj),'transistor')
        NameplateParams.Data = cell(3,3);
    end
    if isequal(class(obj),'inductor')
        NameplateParams.Data = cell(2,3);
    end
    if isequal(class(obj),'capacitor')
        NameplateParams.Data = cell(2,3);
    end

    for i = 1:length(obj.parameters)
        Characteristics.Data{i,1} = obj.parameters(i).name;
        Characteristics.Data{i,2} = char(obj.parameters(i).conditions(:));
        Characteristics.Data{i,3} = obj.parameters(i).min;
        Characteristics.Data{i,4} = obj.parameters(i).typ;
        Characteristics.Data{i,5} = obj.parameters(i).max;
        Characteristics.Data{i,6} = [obj.parameters(i).unit{1}, obj.parameters(i).unit{2}];


        % Create Characteristics
        if isequal(class(obj),'transistor')

            if  strcmp(Characteristics.Data{i,1}, 'Rds')
                if~isempty(Characteristics.Data{i,5})
                    NameplateParams.Data{2,1} = 'Rds(on) - Drain-Source On Resistance';
                    NameplateParams.Data{2,2} = Characteristics.Data{i,5};
                    NameplateParams.Data{2,3} = Characteristics.Data{i,6};
                elseif ~isempty(Characteristics.Data{i,4})
                    NameplateParams.Data{2,1} = 'Rds(on) - Drain-Source On Resistance';
                    NameplateParams.Data{2,2} = Characteristics.Data{i,4};
                    NameplateParams.Data{2,3} = Characteristics.Data{i,6};
                end
            end

            if  strcmp(Characteristics.Data{i,1}, 'Vds')
                if~isempty(Characteristics.Data{i,5})
                    NameplateParams.Data{1,1} = 'Vds - Drain-Source Voltage';
                    NameplateParams.Data{1,2} = Characteristics.Data{i,5};
                    NameplateParams.Data{1,3} = Characteristics.Data{i,6};
                elseif ~isempty(Characteristics.Data{i,4})
                    NameplateParams.Data{1,1} = 'Vds - Drain-Source Voltage';
                    NameplateParams.Data{1,2} = Characteristics.Data{i,4};
                    NameplateParams.Data{1,3} = Characteristics.Data{i,6};
                end
            end

            if  strcmp(Characteristics.Data{i,1}, 'Ids')
                if~isempty(Characteristics.Data{i,5})
                    NameplateParams.Data{3,1} = 'Ids - Drain-Source Current';
                    NameplateParams.Data{3,2} = Characteristics.Data{i,5};
                    NameplateParams.Data{3,3} = Characteristics.Data{i,6};
                elseif ~isempty(Characteristics.Data{i,4})
                    NameplateParams.Data{3,1} = 'Ids - Drain-Source Current';
                    NameplateParams.Data{3,2} = Characteristics.Data{i,4};
                    NameplateParams.Data{3,3} = Characteristics.Data{i,6};
                end
            end

        end


        if isequal(class(obj),'inductor')

            if  strcmp(Characteristics.Data{i,1}, 'L')
                if~isempty(Characteristics.Data{i,5})
                    NameplateParams.Data{1,1} = 'L - Inductance';
                    NameplateParams.Data{1,2} = Characteristics.Data{i,5};
                    NameplateParams.Data{1,3} = Characteristics.Data{i,6};
                elseif ~isempty(Characteristics.Data{i,4})
                    NameplateParams.Data{1,1} = 'L - Inductance';
                    NameplateParams.Data{1,2} = Characteristics.Data{i,4};
                    NameplateParams.Data{1,3} = Characteristics.Data{i,6};
                end
            end

            if  strcmp(Characteristics.Data{i,1}, 'Rdc')
                if~isempty(Characteristics.Data{i,5})
                    NameplateParams.Data{2,1} = 'Rdc - DC Resistance';
                    NameplateParams.Data{2,2} = Characteristics.Data{i,5};
                    NameplateParams.Data{2,3} = Characteristics.Data{i,6};
                elseif ~isempty(Characteristics.Data{i,4})
                    NameplateParams.Data{2,1} = 'Rdc - DC Resistance';
                    NameplateParams.Data{2,2} = Characteristics.Data{i,4};
                    NameplateParams.Data{2,3} = Characteristics.Data{i,6};
                end
            end
        end


        if isequal(class(obj),'capacitor')

            if  strcmp(Characteristics.Data{i,1}, 'C')
                if~isempty(Characteristics.Data{i,5})
                    NameplateParams.Data{1,1} = 'C - Capacitance';
                    NameplateParams.Data{1,2} = Characteristics.Data{i,5};
                    NameplateParams.Data{1,3} = Characteristics.Data{i,6};
                elseif ~isempty(Characteristics.Data{i,4})
                    NameplateParams.Data{1,1} = 'C - Capacitance';
                    NameplateParams.Data{1,2} = Characteristics.Data{i,4};
                    NameplateParams.Data{1,3} = Characteristics.Data{i,6};
                end
            end

            if  strcmp(Characteristics.Data{i,1}, 'Vdc')
                if~isempty(Characteristics.Data{i,5})
                    NameplateParams.Data{2,1} = 'Vdc - DC Rated Voltage';
                    NameplateParams.Data{2,2} = Characteristics.Data{i,5};
                    NameplateParams.Data{2,3} = Characteristics.Data{i,6};
                elseif ~isempty(Characteristics.Data{i,4})
                    NameplateParams.Data{2,1} = 'Vdc - DC Rated Voltage';
                    NameplateParams.Data{2,2} = Characteristics.Data{i,4};
                    NameplateParams.Data{2,3} = Characteristics.Data{i,6};
                end
            end
        end

    end







    %{
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

    %}
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

catch ME
    J=546541698;
end

end

