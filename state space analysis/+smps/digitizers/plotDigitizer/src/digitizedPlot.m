classdef digitizedPlot < handle
    %DIGITIZEDPLOT data storage for digitized image plot
    %   Detailed explanation goes here
    
    properties
        plotData    % x-y data extracted from plot
        logAxes     % [x y], log scale boolean
        normAxes    % [x y], normalized boolean
        
        axisLabels  % string axislabel
        dataLabels  % labels for individual traces
        testConditions ={} %{trace, signal, value, mult, unit}
        SIUnits   
    end
    
    properties (Hidden)
        callingApp
               
        imgData     % in pixels from upper-left of plot area
        normData    % 0-1, in distance from origin   
        
        traceColors % k-means clustering colors identified
        
        pdfPage     % page of the datasheet that plot is on
        pageLoc     % location of the plot on the page
        
        axes        % [Xo Xm Yo Ym] -- limits of the plot
        imgAxes     % [colLeft colRight rowBottom rowTop] -- location of plot on image
        
        img         % RGB data of full plot image
        
        componentType = transistor();

        debug = 0;
    end
    

       
    methods
        function obj = digitizedPlot(app, imgData, pdfLoc, pdfPage)
            %obj = digitizedPlot(app) Construct an instance of this class
            %   Detailed explanation goes here
           
            obj.callingApp = app;
            obj.debug = app.debug;
            if nargin >= 2
                obj.img = imgData;
            end
            if nargin >= 3
                obj.pageLoc = pdfLoc;
            end
            if nargin >= 4
                obj.pdfPage = pdfPage;
            end
        end
        
        function convert(obj, dir)
            %convert(obj, dir) Summary of this method goes here
            %   Detailed explanation goes here
            if nargin == 1
                dir = -1;
            end
            if (dir==1) || (~isempty(obj.imgData) && ~isempty(obj.imgAxes) &&...
                    ~isempty(obj.axes) && isempty(obj.plotData))
               %% Image data and axes calibration data given; create plot data
                for i = 1:length(obj.imgData)
                    ximgdist = obj.imgData{i}(:,1) - obj.imgAxes.colLeft;
                    ximgspan = obj.imgAxes.colRight -  obj.imgAxes.colLeft;
%                     obj.normData{i}(:,1) = 

                    yimgdist = obj.imgAxes.rowBottom - obj.imgData{i}(:,2);
                    yimgspan = obj.imgAxes.rowBottom -  obj.imgAxes.rowTop;
                    
                    obj.normData{i} = [ximgdist/ximgspan, yimgdist/yimgspan];
    
                    %% Convert data, accounting for log Axes
                    if ~(obj.logAxes(1))
                        xspan = obj.axes.Xm - obj.axes.Xo;
                        xPlotData = obj.axes.Xo + obj.normData{i}(:,1)*xspan;
                    else
                        xLogspan = log10(obj.axes.Xm) - log10(obj.axes.Xo);
                        xPlotData = 10.^(log10(obj.axes.Xo) + obj.normData{i}(:,1)*xLogspan);
                    end
                    
                    if ~(obj.logAxes(2))
                        yspan = obj.axes.Ym - obj.axes.Yo;
                        yPlotData = obj.axes.Yo + obj.normData{i}(:,2)*yspan;
                    else
                        yLogspan = log10(obj.axes.Ym) - log10(obj.axes.Yo);
                        yPlotData = 10.^(log10(obj.axes.Yo) + obj.normData{i}(:,2)*yLogspan);
                    end
                    
                    obj.plotData{i} = [xPlotData, yPlotData];
                    
                    
                end
               
            elseif (dir == 2) || (~isempty(obj.plotData) && ~isempty(obj.imgAxes) && ...
                    ~isempty(obj.axes) && isempty(obj.imgData))
                %% Plot data (e.g. grabit) and axes calibration data given; create image data
                for i = 1:length(obj.plotData)
                    if ~(obj.logAxes(1))
                        xdist = obj.plotData{i}(:,1) - obj.axes.Xo;
                        xspan = obj.axes.Xm - obj.axes.Xo;
                        obj.normData{i}(:,1) = xdist/xspan;
                    else
                        xLogspan = log10(obj.axes.Xm) - log10(obj.axes.Xo);
                        obj.normData{i}(:,1) = (log10(obj.plotData{i}(:,1)) - log10(obj.axes.Xo))/xLogspan;
                    end
                    
                    if ~(obj.logAxes(2))
                        ydist = obj.plotData{i}(:,2) - obj.axes.Yo;
                        yspan = obj.axes.Ym - obj.axes.Yo;
                        obj.normData{i}(:,2) = ydist/yspan;
                    else
                        yLogspan = log10(obj.axes.Ym) - log10(obj.axes.Yo);
                        obj.normData{i}(:,2) = (log10(obj.plotData{i}(:,2)) - log10(obj.axes.Yo))/yLogspan;
                    end

                    obj.imgData{i}(:,1) = obj.normData{i}(:,1)* ...
                        (obj.imgAxes.colRight - obj.imgAxes.colLeft) + obj.imgAxes.colLeft;
                    obj.imgData{i}(:,2) = (1-obj.normData{i}(:,2))* ...
                        (obj.imgAxes.rowBottom - obj.imgAxes.rowTop) + obj.imgAxes.rowTop;
                end
            else
                warning('Insufficient data stored to do any conversion');
                return
            end
        end
        
        function loadGrabitData(obj, handle)
            %loadGrabitData(obj, handle) Summary of this method goes here
            %   Detailed explanation goes here
            obj.img = handle.I;
            
            obj.axes.Xo = handle.CalibVals.Xo;
            obj.axes.Xm = handle.CalibVals.Xm;
            obj.axes.Yo = handle.CalibVals.Yo;
            obj.axes.Ym = handle.CalibVals.Ym;
            
            obj.imgAxes.colLeft = handle.CalibVals.Xxo;
            obj.imgAxes.colRight = handle.CalibVals.Xxo + ...
                handle.CalibVals.e1(1)*(handle.CalibVals.Xm-handle.CalibVals.Xo);
            obj.imgAxes.rowBottom = handle.CalibVals.Yyo;
            obj.imgAxes.rowTop = handle.CalibVals.Yyo + ...
                handle.CalibVals.e2(2)*(handle.CalibVals.Ym-handle.CalibVals.Yo);
            
            traces = fieldnames(handle.savedVars);
            for i = 1:length(traces)
                data = eval(['handle.savedVars.' traces{i}]);
                obj.plotData{i} = data;      

%                 obj.imgData{i} = handle.ImDat;
            end
            obj.dataLabels = traces;
            obj.loadAppData();
            obj.convert();

            cOrder = [ 0    0.4470    0.7410;    0.8500    0.3250    0.0980;    0.9290    0.6940    0.1250;    0.4940    0.1840    0.5560;    0.4660    0.6740    0.1880;    0.3010    0.7450    0.9330;    0.6350    0.0780    0.1840];
            obj.traceColors = 255*cOrder(1:size(traces,2), :);

            obj.callingApp.plotImage.Visible = 0;
        end
        
        function loadAppData(obj)
            %loadAppData(obj) Summary of this method goes here
            %   Detailed explanation goes here
            obj.logAxes = [obj.callingApp.LogarithmicXaxisCheckBox.Value, obj.callingApp.LogarithmicYaxisCheckBox.Value];
            obj.normAxes = [obj.callingApp.NormalizedXaxisCheckBox.Value, obj.callingApp.NormalizedYaxisCheckBox.Value];
%             obj.axisLabels = {obj.callingApp.XaxisDropDown.Value, obj.callingApp.YaxisDropDown.Value};
            if ~isempty(obj.callingApp.XaxisTree.SelectedNodes)
                obj.axisLabels{1} = obj.callingApp.XaxisTree.SelectedNodes.Text;
            else
                obj.axisLabels{1} = '';
            end
            if ~isempty(obj.callingApp.YaxisTree.SelectedNodes)
                obj.axisLabels{2} = obj.callingApp.YaxisTree.SelectedNodes.Text;
            else
                obj.axisLabels{2} = '';
            end
            obj.SIUnits = {obj.callingApp.XSIprefixDropDown.Value, obj.callingApp.YSIprefixDropDown.Value};
            
            obj.pdfPage = obj.callingApp.pdfPage();
%             obj.pageLoc = obj.callingApp.plotLoc;
%             if isempty(obj.img)
%                 obj.img = obj.callingApp.plotOrigImage;
%             end
            
%             if isempty(obj.axes)
            	obj.axes.Xo = obj.callingApp.MinXValueEditField.Value;
                obj.axes.Xm = obj.callingApp.MaxXValueEditField.Value;
                obj.axes.Yo = obj.callingApp.MinYValueEditField.Value;
                obj.axes.Ym = obj.callingApp.MaxYValueEditField.Value;
%             end
            
            if isempty(obj.imgAxes)
%                 if length(obj.callingApp.imgAxes) == 4
%                     obj.imgAxes.colLeft = obj.callingApp.imgAxes(1);
%                     obj.imgAxes.colRight = obj.callingApp.imgAxes(2);
%                     obj.imgAxes.rowBottom = obj.callingApp.imgAxes(3);
%                     obj.imgAxes.rowTop = obj.callingApp.imgAxes(4);
%                 else
                    obj.imgAxes.colLeft = [];
                    obj.imgAxes.colRight = [];
                    obj.imgAxes.rowBottom = [];
                    obj.imgAxes.rowTop = [];
%                 end
            end
            
            if isempty(obj.dataLabels) || length(obj.dataLabels) ~= length(obj.imgData)
                if ~isempty(obj.imgData)
                    len = length(obj.imgData);
                elseif ~isempty(obj.plotData)
                    len = length(obj.plotData);
                else
                    return
                end
                obj.dataLabels = cellstr([repmat('Y',len,1), num2str([1:len]')]);
            end
        end
        
        function loadAutoColorDigitizedData(obj, ximg, traces, colors)
            %loadAutoColorDigitizedData(obj, ximg, traces) Summary of this method goes here
            %   Detailed explanation goes here
            obj.loadAppData();
            for i = 1:size(traces,2)
                x = ximg + obj.imgAxes.colLeft;
                y = traces(:,i) + obj.imgAxes.rowTop;
                
                obj.imgData{i} = [x(~isnan(y)), y(~isnan(y))];
            end
            if nargin ==4
                obj.traceColors = colors;
            else
                cOrder = [ 0    0.4470    0.7410;    0.8500    0.3250    0.0980;    0.9290    0.6940    0.1250;    0.4940    0.1840    0.5560;    0.4660    0.6740    0.1880;    0.3010    0.7450    0.9330;    0.6350    0.0780    0.1840];
                obj.traceColors = 255*cOrder(1:size(traces,2), :);
            end
            obj.loadAppData();
            obj.convert(1);
        end
        
        function plotDigitizedData(obj, axes)
            %plotDigitizedData(obj, axes) Summary of this method goes here
            %   Detailed explanation goes here
            hold(axes, 'off');
            imshow(obj.img, 'Parent', axes);
            hold(axes, 'on');
            
            for i = 1:length(obj.imgData)
                if size(obj.traceColors,1) >= i
                    traceColor = obj.traceColors(i,:)/255;
                else 
                    colorWheel = jet;
                    rowInd = round(i/length(obj.imgData)*size(colorWheel,1));
                    traceColor = colorWheel(rowInd,:);
                end

                lines(i) = plot(axes, obj.imgData{i}(:,1), obj.imgData{i}(:,2), ':', 'LineWidth', 3, 'Color', traceColor); 
            end
            
            plot(axes, [obj.imgAxes.colLeft obj.imgAxes.colLeft], [obj.imgAxes.rowBottom obj.imgAxes.rowTop],'-r', 'LineWidth', 3);
            plot(axes, [obj.imgAxes.colRight obj.imgAxes.colRight], [obj.imgAxes.rowBottom obj.imgAxes.rowTop],'-r', 'LineWidth', 3);
            plot(axes, [obj.imgAxes.colLeft obj.imgAxes.colRight], [obj.imgAxes.rowBottom obj.imgAxes.rowBottom],'-r', 'LineWidth', 3);
            plot(axes, [obj.imgAxes.colLeft obj.imgAxes.colRight], [obj.imgAxes.rowTop obj.imgAxes.rowTop],'-r', 'LineWidth', 3);
            legend(axes, lines, obj.dataLabels);
            
            axes.XLim = [0 size(obj.img,2)];
            axes.YLim = [0 size(obj.img,1)];
        end
        
        function addTestCondition(obj, trace, signal, value, mult, unit)
           if isempty(obj.testConditions)
               obj.testConditions = {trace, signal, value, mult, unit};
           else
               obj.testConditions = [obj.testConditions; {trace, signal, value, mult, unit}];
           end
        end
        
        function [str, trace] = getTestConditionStrings(obj, ind)
            if nargin == 1
                start = 1;
                stop = size(obj.testConditions,1);
            elseif nargin == 2
                start = ind;
                stop = ind;
            end
            
            if stop == 0
                str = {''};
                return
            end
            
            for i = start:stop
                unit = obj.testConditions{i,4};
%                 if strcmp(unit, '1')
%                    unit = '';
%                 end
                str{i} = [obj.testConditions{i,2}, '=', num2str(obj.testConditions{i,3}), unit, obj.testConditions{i,5}];
                trace{i} = obj.testConditions{i,1};
            end
            
            if nargin == 2
                str = str{i}; trace = trace{i};
            end
        end
        
        function title = getPlotTitle(obj, len)
            if nargin==1
                len = 'short';
            else
                assert(strcmp(len,'short') || strcmp(len,'long'), 'Valid inputs for input len are ''short'' or ''long''');
            end
            
            if strcmp(len,'short')
                xind = min([length(obj.axisLabels{1}), strfind(obj.axisLabels{1}, ',')-1]);
                yind = min([length(obj.axisLabels{2}), strfind(obj.axisLabels{2}, ',')-1]);
                xName = obj.axisLabels{1}(1:xind);
                yName = obj.axisLabels{2}(1:yind);
            else
                xName = obj.axisLabels{1};
                yName = obj.axisLabels{2};
            end
            
            if obj.logAxes(1)
                xName = ['log ' xName];
            end
            if obj.logAxes(2)
                yName = ['log ' yName];
            end
            
            if obj.normAxes(1)
                xName = ['normalized ' xName];
            end
            if obj.normAxes(2)
                yName = ['normalized ' yName];
            end
            
            if ~isempty(obj.dataLabels)
                vars = {};
                for i = 1:length(obj.dataLabels)
                    ind = min([length(obj.dataLabels{i}), strfind(obj.dataLabels{i}, '=')-1]);
                    vars = {vars{:}, obj.dataLabels{i}(1:ind)};
                end
                zName = [' (' strjoin(unique(vars), ', ') ')'];
            else
                zName = '';
            end
            
            title = [yName '-vs-' xName zName];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        %%          Image Parsing Functions
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function imgAxes = findPlotArea(obj, settings)
            obj.loadAppData();
            
            if nargin < 2
                settings.axisSearchBounds = [0.025 0.2];
            end
            
            bds = settings.axisSearchBounds;
            
            %% Check to see if there is a box around the plot
            bb = ceil(findPlotOuterBox(obj));
            if ~isempty(bb)
                obj.img = obj.img(bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1,:);
            else
            end
            
            % Look for darkest rows/columns of image
            colDark = sum(obj.img, [1,3]);
            rowDark = sum(obj.img, [2,3]);
            
             % Block out outer 2.5% and inner 60%
            colDark(ceil(size(obj.img,2)*bds(2)):floor(size(obj.img,2)*(1-bds(2)))) = inf;
            colDark(1:floor(size(obj.img,2)*bds(1))) = inf;
            colDark(ceil(size(obj.img,2)*(1-bds(1))):size(obj.img,2)) = inf;
            rowDark(ceil(size(obj.img,1)*bds(2)):floor(size(obj.img,1)*(1-bds(2)))) = inf;
            rowDark(1:floor(size(obj.img,1)*bds(1))) = inf;
            rowDark(ceil(size(obj.img,1)*(1-bds(1))):size(obj.img,1)) = inf;
            

            
            % Find rows and cols that look most like horizontal/vertical
            % black lines
            [colc,coli] = mink(colDark, 10);
            [rowc,rowi] = mink(rowDark, 10);
            
            colLeft = coli(coli < size(obj.img,2)/2)';
            if ~isempty(colLeft)
                relColor = colc(coli < size(obj.img,2)/2)/mean(colc(coli < size(obj.img,2)/2));
                colLeft = colLeft(relColor - min(relColor) < 0.05);
                obj.imgAxes.colLeft = min(colLeft);
            else
                obj.imgAxes.colLeft = .25*size(obj.img,2);
            end
            
            colRight = coli(coli > size(obj.img,2)/2)';
            if ~isempty(colRight)
                relColor = colc(coli > size(obj.img,2)/2)/mean(colc(coli > size(obj.img,2)/2));
                colRight = colRight(relColor - min(relColor) < 0.05);
                obj.imgAxes.colRight = max(colRight);
            else
                obj.imgAxes.colRight = .75*size(obj.img,2);
            end
            
            rowTop = rowi(rowi < size(obj.img,1)/2)';
            if ~isempty(colLeft)
                relColor = rowc(rowi < size(obj.img,1)/2)/mean(rowc(rowi < size(obj.img,1)/2));
                rowTop = rowTop(relColor - min(relColor) < 0.05);
                obj.imgAxes.rowTop = min(rowTop);
            else
                obj.imgAxes.rowTop = .25*size(obj.img,1);
            end
            
            rowBottom = rowi(rowi > size(obj.img,2)/2)';
            if ~isempty(rowBottom)
                relColor = rowc(rowi > size(obj.img,1)/2)/mean(rowc(rowi > size(obj.img,1)/2));
                rowBottom = rowBottom(relColor - min(relColor) < 0.05);
                obj.imgAxes.rowBottom = max(rowBottom);
            else
                obj.imgAxes.rowBottom = .75*size(obj.img,1);
            end
            
            imgAxes = [obj.imgAxes.colLeft, obj.imgAxes.colRight, obj.imgAxes.rowBottom, obj.imgAxes.rowTop];
        end
        
        function maskedImg = axesMask(obj)
            try
                locs = [obj.imgAxes.colLeft, obj.imgAxes.colRight, obj.imgAxes.rowBottom, obj.imgAxes.rowTop];
                maskedImg = obj.img;
                if length(locs) == 4
                    locs = round(locs);
                    maskedImg(locs(4)-2:locs(4)+2,locs(1):locs(2),1) = 255;
                    maskedImg(locs(3)-2:locs(3)+2,locs(1):locs(2),1) = 255;
                    maskedImg(locs(4):locs(3),locs(1)-2:locs(1)+2,1) = 255;
                    maskedImg(locs(4):locs(3),locs(2)-2:locs(2)+2,1) = 255;

                    maskedImg(locs(4)-2:locs(4)+2,locs(1):locs(2),2:3) = 0;
                    maskedImg(locs(3)-2:locs(3)+2,locs(1):locs(2),2:3) = 0;
                    maskedImg(locs(4):locs(3),locs(1)-2:locs(1)+2,2:3) = 0;
                    maskedImg(locs(4):locs(3),locs(2)-2:locs(2)+2,2:3) = 0;
                end
            catch
                x=1;
            end
            
        end
        
        function maskedImg = OCRaxes(obj, mask, settings)
            if nargin < 3
                settings.labelSearchRadius = round(min(size(obj.img,1), size(obj.img,2))*.1);
            end
            if nargin<2
                mask = obj.img;
            end
            r = settings.labelSearchRadius;
            
            obj.resetPlotLabels();

            % Regions for axis labels
            rois = [1, obj.imgAxes.rowBottom, size(obj.img,2)-2, size(obj.img,1)-obj.imgAxes.rowBottom-1;  %
                1, 1, obj.imgAxes.colLeft, size(obj.img,1)-2;
                ];

            % Check for recognizeable text
            ocrResults = ocr(obj.img, rois);
            settings.ocrDilation = 0;
%             newocrResults = obj.tunedOCR(rois, obj.img, 0, settings);

            % For Y-axis label, rotate the image to make text horizontal
            ocrResults = [ocrResults(1); ...
                ocr(rot90(obj.img,3), [rois(2,2) rois(2,1) rois(2,4), rois(2,3)])];
%             newocrResults = [newocrResults(1); ...
%                 obj.tunedOCR([rois(2,2) rois(2,1) rois(2,4), rois(2,3)], rot90(obj.img,3), 0, settings)];

            
%             ocrResults = newocrResults;
            
            %% Check for identified axis labels
%             Xlabel = intersect(ocrResults(1).Words, obj.callingApp.XaxisDropDown.Items);
%             Ylabel = intersect(ocrResults(2).Words, obj.callingApp.YaxisDropDown.Items);
            possibleYLabels = regexprep(ocrResults(2).Words, '[.,''!?-><]', '');
            possibleXLabels = regexprep(ocrResults(1).Words, '[.,''!?-><]', '');
            [Xlabel, ~, XLloc] = intersect(possibleXLabels, {obj.callingApp.XaxisTree.Children.Text});
            [Ylabel, ~, YLloc] = intersect(possibleYLabels, {obj.callingApp.YaxisTree.Children.Text});

            
            
            if ~isempty(Xlabel)
%                 obj.callingApp.XaxisDropDown.Value = Xlabel;
%                 obj.callingApp.XaxisDropDown.FontColor = [1 0 0];
                obj.callingApp.XaxisTree.SelectedNodes = obj.callingApp.XaxisTree.Children(XLloc);
                obj.callingApp.XaxisTree.Children(XLloc).expand();
                obj.callingApp.XaxisTree.scroll(obj.callingApp.XaxisTree.Children(XLloc));
                obj.callingApp.XaxisTree.FontColor = [0 .7 0];
            end

            if ~isempty(Ylabel)
%                 obj.callingApp.YaxisDropDown.Value = Ylabel;
%                 obj.callingApp.YaxisDropDown.FontColor = [1 0 0];
                obj.callingApp.YaxisTree.SelectedNodes = obj.callingApp.YaxisTree.Children(YLloc);
                obj.callingApp.YaxisTree.Children(YLloc).expand();
                obj.callingApp.YaxisTree.scroll(obj.callingApp.YaxisTree.Children(YLloc));
                obj.callingApp.YaxisTree.FontColor = [0 .7 0];
            end
            
            %% Check for identified unit prefixes
            Xunit = ocrResults(1).Words(~cellfun(@isempty,regexp(ocrResults(1).Words,'\(.+\)')));
            if ~isempty(Xunit)
                Xunit = strrep(Xunit, '(','');
                Xunit = strrep(Xunit, ')','');
                Xunit = Xunit{end}(1);
                Xunit = intersect(obj.callingApp.XSIprefixDropDown.Items, Xunit);
            end

            Yunit = ocrResults(2).Words(~cellfun(@isempty,regexp(ocrResults(2).Words,'\(.+\)')));
            if ~isempty(Yunit)
                Yunit = strrep(Yunit, '(','');
                Yunit = strrep(Yunit, ')','');
                Yunit = Yunit{end}(1);
                Yunit = intersect(obj.callingApp.XSIprefixDropDown.Items, Yunit);
            end

            if ~isempty(Xunit)
                obj.callingApp.XSIprefixDropDown.Value = Xunit;
                obj.callingApp.XSIprefixDropDown.FontColor = [1 0 0];
            end

            if ~isempty(Yunit)
                obj.callingApp.YSIprefixDropDown.Value = Yunit;
                obj.callingApp.YSIprefixDropDown.FontColor = [1 0 0];
            end
            
            %% Look for axis values
            % Approx locations of minimum/maximum axis labels
            % [x y width height]
            rois = [obj.imgAxes.colLeft-r-3, obj.imgAxes.rowTop-r/2, r, r;
                obj.imgAxes.colLeft-r-3, obj.imgAxes.rowBottom-r/2, r, r;
                obj.imgAxes.colLeft-r/2, obj.imgAxes.rowBottom+3, r, r;
                obj.imgAxes.colRight-r/2, obj.imgAxes.rowBottom+3, r, r;
                ];
            
            % Keep everything within the image extents
            len = size(rois,1);
            rois = min(rois, [size(obj.img,2)*ones(len,1) - rois(:,3), ...
                size(obj.img,1)*ones(len,1) - rois(:,4), inf*ones(len,2) ]);
            rois = max(rois, ones(len,4) );
            
            % Complete axis
            rois = [rois;
                obj.imgAxes.colLeft-3, obj.imgAxes.rowBottom+3, size(obj.img,2)-obj.imgAxes.colLeft, size(obj.img,1)-(obj.imgAxes.rowBottom+3);  %
                1, 1, obj.imgAxes.colLeft-5, size(obj.img,1)-2;
                ];
            
            %% If we did find the axis label, adjust the boxes for numerical ticks
            if ~isempty(Xlabel)
                xlabelLoc = ocrResults(1).WordBoundingBoxes(strcmp(possibleXLabels, Xlabel),:);
                xlabelEnd = xlabelLoc(2)-5;
                rois(5,4) = xlabelEnd-rois(5,2);
            end

            if ~isempty(Ylabel)
                ylabelLoc = ocrResults(2).WordBoundingBoxes(strcmp(possibleYLabels, Ylabel),:);
                ylabelEnd = ylabelLoc(2) + ylabelLoc(4)+5;
                rois(6,:) = [ylabelEnd, 1, rois(6,3)-ylabelEnd, rois(6,4)];
            end
            
            % Check for recognizeable text
            ocrResults = ocr(obj.img, rois,'CharacterSet','.0123456789');
            
            w = warning ('off','all');
                for i = 1:length(ocrResults)
                    structocrResults(i) = struct(ocrResults(i));
                end
                ocrResults = structocrResults;
            warning(w)
            
            if(0)
                roiP = rois(6,:);
                imshow(obj.img(roiP(2):roiP(2)+roiP(4), roiP(1):roiP(1)+roiP(3) ))
                results = tunedOCR(obj, roiP, obj.img, 1)
            end
            
            % If we found nothing, try again with tuned OCR
            for i = 1:length(ocrResults)
                if isempty(ocrResults(i).Text)
                    tOCR = tunedOCR(obj, rois(i,:), obj.img, 1);
                    if ~isempty(tOCR)
                        ocrResults(i) = tOCR;
                    end
                end
            end
            
%             for i = 1:length(ocrResults)
%                 if isempty(ocrResults(i).Words)
%                     newocrResults2(i) = obj.tunedOCR(rois(i,:), obj.img, 1);
%                 else
%                     newocrResults2(i).Words = ocrResults(i).Words;
%                     newocrResults2(i).WordBoundingBoxes = ocrResults(i).WordBoundingBoxes;
%                 end
%             end
%             
%             ocrResults = newocrResults2;


            %% Check for identified axis ticks anywhere
            xaxisvals = cellfun(@(x) str2double(x), ocrResults(5).Words,'UniformOutput', false);
            xaxisvals = [xaxisvals{:}];

            yaxisvals = cellfun(@(x) str2double(x), ocrResults(6).Words,'UniformOutput', false);
            yaxisvals = [yaxisvals{:}];
            

            
            %outlier identification
            yOL = yaxisvals < mean(yaxisvals) - 1.5*iqr(yaxisvals) | yaxisvals > mean(yaxisvals) + 1.5*iqr(yaxisvals);
            xOL = xaxisvals < mean(xaxisvals) - 1.5*iqr(xaxisvals) | xaxisvals > mean(xaxisvals) + 1.5*iqr(xaxisvals);
            
            % check for log axis
            if length(yaxisvals) >= 3
                if std(diff(log10(yaxisvals(~yOL))))*5 < std(diff(yaxisvals(~yOL)))
                    obj.callingApp.LogarithmicYaxisCheckBox.Value = 1;
                    logY = log10(yaxisvals);
                    yOL = logY < mean(logY) - 1.5*iqr(logY) | logY > mean(logY) + 1.5*iqr(logY);
                end
            end
            
            if length(xaxisvals) >= 3
                if std(diff(log10(xaxisvals(~xOL))))*5 < std(diff(xaxisvals(~xOL)))
                    obj.callingApp.LogarithmicXaxisCheckBox.Value = 1;
                    logX = log10(xaxisvals);
                    xOL = logX < mean(logX) - 1.5*iqr(logX) | logX > mean(logX) + 1.5*iqr(logX);
                end
            end
            
            xaxisvals(xOL) = [];
            yaxisvals(yOL) = [];

            alternates = [max(yaxisvals), min(yaxisvals), min(xaxisvals), max(xaxisvals)];
            
            %% Try at looking at all labels and doing a linear regression
            try
                xaxislocs = [ocrResults(5).WordBoundingBoxes(:,1) + ocrResults(5).WordBoundingBoxes(:,3)/2, ...
                    ocrResults(5).WordBoundingBoxes(:,2) + ocrResults(5).WordBoundingBoxes(:,4)/2];
                xaxislocs(xOL,:) = [];
                [P,S, mu] = polyfit(xaxislocs(~isnan(xaxisvals),1),xaxisvals(~isnan(xaxisvals))',1);
%                 fiterr =  xaxisvals(~isnan(xaxisvals)) - polyval(P,xaxislocs(~isnan(xaxisvals),1),S, mu);
                xlims = polyval(P,[obj.imgAxes.colLeft, obj.imgAxes.colRight],S, mu);
                xlims = round(xlims,3);
            catch
                if(debug)
                    warning('x-axis label regression failed');
                end
            end
            try
                yaxislocs = [ocrResults(6).WordBoundingBoxes(:,1) + ocrResults(6).WordBoundingBoxes(:,3)/2, ...
                    ocrResults(6).WordBoundingBoxes(:,2) + ocrResults(6).WordBoundingBoxes(:,4)/2];
                yaxislocs(yOL,:) = [];
                ws = warning('off','all');  % Turn off warning
                if ~obj.callingApp.LogarithmicYaxisCheckBox.Value
                    [P,S, mu] = polyfit(yaxislocs(~isnan(yaxisvals),2),yaxisvals(~isnan(yaxisvals))',1);
                    ylims = polyval(P,[obj.imgAxes.rowTop, obj.imgAxes.rowBottom],S, mu);
                    ylims = round(ylims,3);
                else
                    [P,S, mu] = polyfit(yaxislocs(~isnan(yaxisvals),2),log10(yaxisvals(~isnan(yaxisvals)))',1);
                    ylims = polyval(P,[obj.imgAxes.rowTop, obj.imgAxes.rowBottom],S, mu);
                    ylims = 10.^round(ylims,3);
                end
                warning(ws);
            catch
                if(debug)
                    warning('y-axis label regression failed');
                end
            end

            %% check if we found axis labels at the ends of the axis
            useable = cellfun(@(x) ~isempty(x), {ocrResults(1:4).Words});
            % also check if it is consistent with tick labels
            if length(xaxisvals) >= 3
                if ~isempty(ocrResults(3).Words)
                    useable(3) = str2double(ocrResults(3).Words{1}) > mean(xaxisvals) - 1.5*iqr(xaxisvals) & str2double(ocrResults(3).Words{1}) > mean(xaxisvals) + 1.5*iqr(xaxisvals);
                    if exist('xlims','var')
                        useable(3) = useable(3) && abs(str2double(ocrResults(3).Words{1}) - xlims(1)) < iqr(xaxisvals);
                    end
                end
                if ~isempty(ocrResults(4).Words)
                    useable(4) = str2double(ocrResults(4).Words{1}) > mean(xaxisvals) - 1.5*iqr(xaxisvals) & str2double(ocrResults(4).Words{1}) < mean(xaxisvals) + 1.5*iqr(xaxisvals);
                    if exist('xlims','var')
                        useable(4) = useable(4) && abs(str2double(ocrResults(4).Words{1}) - xlims(2)) < iqr(xaxisvals);
                    end
                end
            end
            
            if length(yaxisvals) >= 3
                if ~isempty(ocrResults(1).Words{1})
                    if ~obj.callingApp.LogarithmicYaxisCheckBox.Value
                        useable(1) = str2double(ocrResults(1).Words{1}) > mean(yaxisvals) - 1.5*iqr(yaxisvals) & str2double(ocrResults(1).Words{1}) < mean(yaxisvals) + 1.5*iqr(yaxisvals);
                    else
                        yLog = log10(yaxisvals);
                        useable(1) = log10(str2double(ocrResults(1).Words{1})) > mean(yLog) - 1.5*iqr(yLog) & log10(str2double(ocrResults(1).Words{1})) < mean(yLog) + 1.5*iqr(yLog);
                    end
                    if exist('ylims','var')
                        useable(1) = useable(1) && abs(str2double(ocrResults(1).Words{1}) - ylims(1)) < iqr(yaxisvals);
                    end
                end
                if ~isempty(ocrResults(2).Words{1})
                    if ~obj.callingApp.LogarithmicYaxisCheckBox.Value
                        useable(2) = str2double(ocrResults(2).Words{1}) > mean(yaxisvals) - 1.5*iqr(yaxisvals) & str2double(ocrResults(2).Words{1}) < mean(yaxisvals) + 1.5*iqr(yaxisvals);
                    else
                        yLog = log10(yaxisvals);
                        useable(2) = log10(str2double(ocrResults(2).Words{1})) > mean(yLog) - 1.5*iqr(yLog) & log10(str2double(ocrResults(2).Words{1})) < mean(yLog) + 1.5*iqr(yLog);
                    end
                    if exist('ylims','var')
                        useable(2) = useable(2) && abs(str2double(ocrResults(2).Words{1}) - ylims(2)) < iqr(yaxisvals);
                    end
                end
            end
            
            %% Decide on which axis limits to use
            for i = 1:4
                if useable(i) == 1
                    % Found value in the right location
                    val = str2double(ocrResults(i).Words{1});
                elseif (i == 1 || i == 2) && exist('ylims','var')
                    % Use linear regression results
                    val = ylims(i);
                elseif (i == 3 || i == 4) && exist('xlims','var')
                    % Use linear regression results
                    val = xlims(i-2);
                elseif length(alternates) == 4
                    % Use the closest number we found
                    val = alternates(i);   
                else
                    val = nan;
                end

                if ~isnan(val)
                    switch i
                        case 1
                            obj.callingApp.MaxYValueEditField.Value = val;
                            obj.callingApp.MaxYValueEditField.FontColor = [1 0 0];
                        case 2
                            obj.callingApp.MinYValueEditField.Value = val;
                            obj.callingApp.MinYValueEditField.FontColor = [1 0 0];
                        case 3
                            obj.callingApp.MinXValueEditField.Value = val;
                            obj.callingApp.MinXValueEditField.FontColor = [1 0 0];
                        case 4
                            obj.callingApp.MaxXValueEditField.Value = val;
                            obj.callingApp.MaxXValueEditField.FontColor = [1 0 0];
                    end
                end

            end


            Iocr = insertObjectAnnotation(mask, 'rectangle', ...
                       rois, {'Ymax', 'Ymin', 'Xmin', 'Xmax', 'XTicks', 'YTicks'}, ...
                        'Color', 'red');
            maskedImg = Iocr;
            if sum(useable) > 0
%                         validResults = ocrResults(useable ==1);
%                         Iocr = insertObjectAnnotation(Iocr, 'rectangle', ...
%                                    validResults.WordBoundingBoxes, ...
%                                    validResults.WordConfidences);
%                         app.plotMask = double(Iocr) - double(app.plotOrigImage);
%                         
            end
        end
        
        function results = tunedOCR(obj, rois, img, isNumeric, settings)
           %https://www.mathworks.com/matlabcentral/answers/377444-why-ocr-function-doesn-t-recognize-the-numbers
            % Localize words
            if nargin < 3
                img = obj.img;
            end
            
            if nargin < 4
                isNumeric = 0;
            end
            
            if nargin < 5
                settings = [];
            end
            
            if ~isfield(settings, 'boundingBoxDilation')
                settings.boundingBoxDilation = 4;
            end
            if ~isfield(settings, 'ocrDilation')
                settings.ocrDilation = 1;
            end
            
            rois = round(rois);
                       
            %% dilate image to improve clarity
            img = 255-imdilate((255-img),strel('disk',settings.ocrDilation));
            
            BW = imbinarize(rgb2gray(img));
            
            
            bboxes = [];
            for i = 1:size(rois,1)
                subimg = BW(rois(i,2):rois(i,2) + rois(i,4),rois(i,1):rois(i,1) + rois(i,3));
                subimg = imdilate(~subimg,strel('disk',settings.boundingBoxDilation));
                s = regionprops(subimg,'BoundingBox');
                if ~isempty(s)
                    bboxes = [bboxes; vertcat(s(:).BoundingBox) + [rois(i,1:2) 0 0]];
                end
            end
            
            if(0)
                % plot bboxes -- Debug section
                BW1 = imdilate(~BW,strel('disk',settings.boundingBoxDilation));
               
                labels = num2cell(1:size(bboxes,1));
                labels = cellstr(num2str([labels{:}]'));
                labeledImg = insertObjectAnnotation(uint8(255*(~BW1)), 'rectangle',bboxes, labels, 'Color', 'red');
                
                labels = num2cell(1:size(rois,1));
                labels = cellstr(num2str([labels{:}]'));
                labeledImg = insertObjectAnnotation(labeledImg, 'rectangle',rois, labels, 'Color', 'blue');
                
                imshow(labeledImg)
            end
            
            
%             % Pre-process image to make letters thicker
%             BW = imdilate(~BW,strel('disk',1,0));
            
            % Sort boxes into rois
            w = warning ('off','all');
                results = struct(ocrText.empty);
            warning(w)
            for i = 1:size(rois,1)
                if ~isempty(bboxes)
                    cont = diag(rectint(rois(i,:), bboxes) > .5*diag(rectint(bboxes,bboxes)));
                    if isNumeric
                        ocrResults = ocr(img,bboxes(cont,:),'CharacterSet','.0123456789','TextLayout','word');
                    else
                        ocrResults = ocr(img,bboxes(cont,:),'TextLayout','word');
                    end
                    results(i).Text = [ocrResults.Text];
                    results(i).Words = [ocrResults(:).Words];
                    results(i).WordBoundingBoxes = cat(1,ocrResults(:).WordBoundingBoxes);
                    results(i).WordConfidences = cat(1,ocrResults(:).WordConfidences);
                    results(i).CharacterBoundingBoxes = cat(1,ocrResults(:).CharacterBoundingBoxes);
                    results(i).CharacterConfidences = cat(1,ocrResults(:).CharacterConfidences);
                end
            end
                       
        end
        
        function bb = findPlotOuterBox(obj, settings)
            if nargin < 2
                settings.pixelAverageWidth = 5;
                settings.pixelAverageThreshold = 0.25;
                settings.lineClusterThreshold = 100;
                settings.imageOuterMargin = 50;
            end

            imgPlot = obj.img;%(obj.imgAxes.rowTop-settings.imageOuterMargin:obj.imgAxes.rowBottom+settings.imageOuterMargin,...
                    %obj.imgAxes.colLeft-settings.imageOuterMargin:obj.imgAxes.colRight+settings.imageOuterMargin,:);
            
            BWdata = imbinarize(rgb2gray(imgPlot));

            stats = regionprops(~BWdata,'FilledArea', 'BoundingBox');
            
            refArea = size(imgPlot,1)*size(imgPlot,2);
            
            areaErr = abs([stats.FilledArea] - refArea)/refArea; 
            
            [minErr, loc] = min(areaErr);
            if minErr < .15
                bb = stats(loc).BoundingBox;
            else
                bb = [];
            end
        end
        
        function bboxes = findSimilarPlots(obj, datasheet)
            
%             img = imbinarize(rgb2gray(obj.img));
            BWdata = imbinarize(rgb2gray(datasheet));
            
            bb = ceil(findPlotOuterBox(obj));
            if isempty(bb)
                refArea = abs(obj.imgAxes.rowTop-obj.imgAxes.rowBottom)*abs(obj.imgAxes.colLeft-obj.imgAxes.colRight);
            else
                refArea = bb(3)*bb(4);
            end
%             stats = regionprops(~BWdata)
            stats = regionprops(~BWdata,'FilledArea', 'BoundingBox');
            
            areaErr = abs([stats.FilledArea] - refArea)/refArea;
            otherPlots = stats(areaErr < 0.1);
            
            if isempty(bb)
                for i = 1:length(otherPlots)
                    otherPlots(i).BoundingBox = [otherPlots(i).BoundingBox] - [obj.imgAxes.colLeft, obj.imgAxes.rowTop, 0 ,0];
                    otherPlots(i).BoundingBox(3:4) = [size(obj.img,2), size(obj.img,1)];
                end
            else
%                 for i = 1:length(otherPlots)
%                     otherPlots(i).BoundingBox = [otherPlots(i).BoundingBox] - [obj.imgAxes.colLeft, obj.imgAxes.rowTop, 0 ,0];
%                     otherPlots(i).BoundingBox(3:4) = [size(obj.img,2), size(obj.img,1)];
%                 end
            end
            
            bboxes = vertcat(otherPlots.BoundingBox);
            
%             labels = cellstr(num2str([1:size(bboxes,1)]'))'
%             Iocr = insertObjectAnnotation(datasheet, 'rectangle', ...
%                         bboxes, labels, ...
%                         'Color', 'red');
% 
%             imshow(Iocr)

            
%             leftBound = img(obj.imgAxes.rowTop:obj.imgAxes.rowBottom,obj.imgAxes.colLeft,:);
%             rightBound = img(obj.imgAxes.rowTop:obj.imgAxes.rowBottom,obj.imgAxes.colRight,:);
%             topBound = img(obj.imgAxes.rowTop,obj.imgAxes.colLeft:obj.imgAxes.colRight,:);
%             bottomBound = img(obj.imgAxes.rowBottom,obj.imgAxes.colLeft:obj.imgAxes.colRight,:);
%             
%             err = nan*ones(size(datasheet,2)-size(topBound,2), size(datasheet,1)-size(leftBound,1));
%             
%             mask = boolean(zeros(size(leftBound,1), size(topBound,2), 1));
%             mask(1,:,:) = 1;  mask(end,:,:) = 1;
%             mask(:,1,:) = 1;  mask(:,end,:) = 1;
%             
%             ref = boolean(zeros(size(leftBound,1), size(topBound,2), 1));
%             ref(:,1,:) = leftBound;  ref(:,end,:) = rightBound;
%             ref(1,:,:) = topBound;  ref(end,:,:) = bottomBound;
%             
%             origLoc = obj.pageLoc(1:2) + [obj.imgAxes.colLeft, obj.imgAxes.rowTop];
%             j = origLoc(2)-1; i = origLoc(1)-1;
%             origSelection = BWdata(j:j+size(leftBound,1)-1,i:i+size(topBound,2)-1,:);
%             refErr = sum(abs(origSelection.*mask - ref), [1:2]);
            
%             h = waitbar(0, 'Checking for other plots');
%             try
%                 for i = 1:4:size(datasheet,2)-size(topBound,2)
%                     for j = 1:4:size(datasheet,1)-size(leftBound,1)
%                         selection = BWdata(j:j+size(leftBound,1)-1,i:i+size(topBound,2)-1,:);
%                         err(i,j) = sum(abs(selection.*mask - ref), [1:2]);
%                     end
%                     waitbar(i/(size(datasheet,2)-size(topBound,2)),h); 
%                 end
%             catch
%                 err;
%             end
%             close(h)
            
        end
        
        function resetPlotLabels(obj)
            obj.callingApp.XaxisDropDown.Value = "<not selected>";
            obj.callingApp.XaxisDropDown.FontColor = [0 0 0];
            obj.callingApp.YaxisDropDown.Value = "<not selected>";
            obj.callingApp.YaxisDropDown.FontColor = [0 0 0];
            
            obj.callingApp.XaxisTree.SelectedNodes = obj.callingApp.XaxisTree.Children(1);
            obj.callingApp.XaxisTree.FontColor = [0 0 0];
            obj.callingApp.XaxisTree.collapse('all');
            obj.callingApp.YaxisTree.SelectedNodes = obj.callingApp.YaxisTree.Children(1);
            obj.callingApp.YaxisTree.FontColor = [0 0 0];
            obj.callingApp.YaxisTree.collapse('all');
            
            obj.callingApp.MaxYValueEditField.Value = 0;
            obj.callingApp.MaxYValueEditField.FontColor = [0 0 0];
            obj.callingApp.MinYValueEditField.Value = 0;
            obj.callingApp.MinYValueEditField.FontColor = [0 0 0];
            obj.callingApp.MinXValueEditField.Value = 0;
            obj.callingApp.MinXValueEditField.FontColor = [0 0 0];
            obj.callingApp.MaxXValueEditField.Value = 0;
            obj.callingApp.MaxXValueEditField.FontColor = [0 0 0];
            obj.callingApp.YSIprefixDropDown.Value = '';
            obj.callingApp.YSIprefixDropDown.FontColor = [0 0 0];
            obj.callingApp.XSIprefixDropDown.Value = '';
            obj.callingApp.XSIprefixDropDown.FontColor = [0 0 0];
            
            obj.callingApp.NormalizedYaxisCheckBox.Value = 0;
            obj.callingApp.NormalizedXaxisCheckBox.Value = 0;
            obj.callingApp.LogarithmicYaxisCheckBox.Value = 0;
            obj.callingApp.LogarithmicXaxisCheckBox.Value = 0;
        end
        
        function nColors = imgTraceColors(obj, imgPlot, pctCoverage)
            if nargin == 1
                pctCoverage = 0.99;
                imgPlot = obj.removeUncoloredPixels();
            elseif nargin == 2
                pctCoverage = 0.99;
            end
            
            %%
            hsv = rgb2hsv(imgPlot);
            [binCount,~] = histcounts(hsv(:,:,1));
            binCount = sort(binCount, 'descend');
            binCoverage = cumsum(binCount)/(size(hsv,1)*size(hsv,2));
            nColors = find(binCoverage > pctCoverage,1,'first');

            nColors = min(max(nColors,2),10);  
            
        end
        
        function imgPlot = removeUncoloredPixels(obj, threshold)
            if nargin == 1
                threshold = 0.05;
            end
            imgPlot = obj.img(obj.imgAxes.rowTop:obj.imgAxes.rowBottom, obj.imgAxes.colLeft:obj.imgAxes.colRight,:);
            hsv = rgb2hsv(imgPlot);
            greyscale = uint8(hsv(:,:,2) < threshold);
            imgPlot(:,:,1) = imgPlot(:,:,1) + (255 - imgPlot(:,:,1)).*greyscale;
            imgPlot(:,:,2) = imgPlot(:,:,2) + (255 - imgPlot(:,:,2)).*greyscale;
            imgPlot(:,:,3) = imgPlot(:,:,3) + (255 - imgPlot(:,:,3)).*greyscale;
        end
        
        function colorClustering(obj, nColors, settings)
            obj.loadAppData();

            imgPlot = obj.removeUncoloredPixels();
            
            if nargin < 2
                nColors = obj.imgTraceColors(imgPlot);
            end
            
            if nargin < 3
                settings.kMeansNAttempts = 3;
                settings.kMeansThreshold = 1e-6;
            end
            
            pixel_labels = imsegkmeans(imgPlot,nColors,'NumAttempts', ...
                settings.kMeansNAttempts, 'Threshold', settings.kMeansThreshold);
            
            colors = zeros(nColors,3);
            timesChanged = zeros(nColors, size(imgPlot,2));
            for i = 1:nColors
                R = imgPlot(:,:,1);
                G = imgPlot(:,:,2);
                B = imgPlot(:,:,3);
                colors(i,:) = [ mean(R(pixel_labels==i)), mean(G(pixel_labels==i)), mean(B(pixel_labels==i))];
                timesChanged(i,:) = sum(abs(diff(pixel_labels == i,1)),1);
            end
            [~,worstCandidates] = sort(sum((timesChanged == 2) - (timesChanged > 2),2));
            bestCandidates = flip(worstCandidates);
            
            white = find(sum(colors,2) > 230*3);
            if ~isempty(white)
                [~, IA, ~] = intersect(bestCandidates, white);
                bestCandidates(IA) = [];
            end
            
            colors = colors(bestCandidates,:);
            
            ximg = round(linspace(1, size(imgPlot,2), 50))';
            trace = nan*ones(50,length(bestCandidates));
            for i = 1:length(bestCandidates)
                src = pixel_labels==bestCandidates(i);
                filt = movmean(src,5);
                
                %% convolution with a triangular pulse 
                % to find middle of largest section of line as start point
                % --helps to avoid noise and legends
                preInt = round(size(imgPlot,2)/2);
                h = 1:preInt;
                h = [h, flip(h)]/max(h);
                
                goodLoc = 1*(timesChanged(bestCandidates(i),:)'==2);
                
                
                % Seems like sometimes this is a row, sometimes a column...
                if (size(goodLoc,1) == 1) && (size(h,1) == 1)
                    filtFind = conv(h,goodLoc);
                    filtFind = filtFind(preInt:end-preInt).*goodLoc;
                else
                    filtFind = conv(h,goodLoc');
                    filtFind = filtFind(preInt:end-preInt).*goodLoc';
                end
                
                % start from wherever we're most confident we've found the
                % line
                startLoc = find(filtFind(ximg) == max(filtFind(ximg)),1);
                
                yimg = nan*ones(1, 50);
                xloc = ximg(startLoc);
                loc = find(filt(:,xloc) == max(filt(:,xloc)));
                if length(loc) > 1
                    loc = loc(round(length(loc)/2));
                    yimg(startLoc) = loc;
                elseif length(loc) == 1
                    yimg(startLoc) = loc;
                elseif isempty(loc) || max(filt(:,xloc)) < 0.01
                    yimg(startLoc) = NaN;
                end
                
                 %% going right
                for j = startLoc+1:49
                   
                    xloc = ximg(j);
                    loc = find((filt(:,xloc) == max(filt(:,xloc))) & (filt(:,xloc)>0));
                    yjump = abs(loc - yimg(j-1));
                    if(j>2)
                        if ~isnan(yimg(j-2))
                            % if data is available, use slope projection
                            yjump = abs( loc - 2*yimg(j-1) + yimg(j-2) );  
                        end
                    end
                    
                    if length(loc) > 1 && ~isnan(yimg(j-1))
                        if min(yjump) < size(imgPlot,1)*.2
                            yimg(j) = loc(yjump == min(yjump));
                        end
                    elseif length(loc) == 1 && ~isnan(yimg(j-1))
                        if min(yjump) < size(imgPlot,1)*.2
                            yimg(j) = loc;
                        end
                    else
                        continue;
                    end
                end
                
                %% going left
                for j = flip(1:startLoc)
                    xloc = ximg(j);
                    loc = find((filt(:,xloc) == max(filt(:,xloc))) & (filt(:,xloc)>0));
                    yjump = abs(loc - yimg(j+1));
                    if j < length(ximg)-2
                        if ~isnan(yimg(j+2))
                            % if data is available, use slope projection
                            yjump = abs( loc - 2*yimg(j+1) + yimg(j+2) ); 
                        end
                    end
                    
                    if length(loc) > 1 && ~isnan(yimg(j+1))
                        if min(yjump) < size(imgPlot,1)*.2
                            yimg(j) = loc(yjump == min(yjump));
                        end
                    elseif length(loc) == 1 && ~isnan(yimg(j+1))
                        if min(yjump) < size(imgPlot,1)*.2
                            yimg(j) = loc;
                        end
                    else
                        continue;
                    end
                end
                trace(:,i) = yimg;
            end
            if ~isempty(trace)
                obj.loadAutoColorDigitizedData(ximg, trace, colors);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        %%          Greyscale Image Parsing Functions
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [nTraces, L] = BWimgTraceNumber(obj, imgPlot, varargin)
            
            settings.pixelAverageWidth = 5;
            settings.pixelAverageThreshold = 0.25;
            settings.lineClusterThreshold = 100;
            settings.imageInnerMargin = 10;
            
            if nargin >= 3
                assert(mod(numel(varargin),2)==0, 'Additional parameters must come in pairs');
                for i = 1:2:numel(varargin)
                    assert(ischar(varargin{i}), 'properties must be character vectors');
                    assert(isnumeric(varargin{i+1}), 'values must be numeric');
                    settings.(varargin{i}) = varargin{i+1};
                end
            end
            
            if nargin == 1 || isempty(imgPlot)
                imgPlot = obj.img(obj.imgAxes.rowTop+settings.imageInnerMargin:obj.imgAxes.rowBottom-settings.imageInnerMargin,...
                    obj.imgAxes.colLeft+settings.imageInnerMargin:obj.imgAxes.colRight-settings.imageInnerMargin,:);
            end
            
            BW = imbinarize(rgb2gray(imgPlot));
            BW = movmean(movmean(BW,settings.pixelAverageWidth),settings.pixelAverageWidth,2)<settings.pixelAverageThreshold ;
            L = bwlabel(BW);
            
            for i = 1:max(max(L))
                count = sum(L==i,[1,2]);
                if count < settings.lineClusterThreshold 
                    L(L==i) = 0;
                end
            end
            
            timesChanged = sum(diff(L~=0,1,1)>0,1);
            [N,edges] = histcounts(timesChanged, 'BinMethod','integer');
            if edges(1) == -0.5
                edges = edges(2:end);
                N = N(2:end);
            end
            maxBinInd = find(N==max(N));
            nTraces = (edges(maxBinInd) + edges(maxBinInd+1))/2;  
            
        end
        
        function BWlabel(obj, nTraces, L, settings)
            if nargin < 4
                settings.imageInnerMargin = 10;
                settings.pixelAverageWidth = 5;
                settings.pixelAverageThreshold = 0.25;
                settings.nXpoints = 50;
                settings.maxYRelPtJump = 0.2;
            end
            
            if nargin == 1
               [nTraces, L] = BWimgTraceNumber(obj);
               nFound = nTraces;
            elseif nargin >= 2
               [nFound, L] = BWimgTraceNumber(obj);
            end 
            
            %% Check if default filtering lost some traces
            [nTracesNF, LNF] = obj.BWimgTraceNumber([], 'pixelAverageWidth',1);
            if nTracesNF > nFound
                nTraces = nTracesNF;
                L = LNF;
                settings.pixelAverageWidth = 1;
            end
            

            imgPlot = obj.img(obj.imgAxes.rowTop:obj.imgAxes.rowBottom, obj.imgAxes.colLeft:obj.imgAxes.colRight,:);
            
            %Re-insert the inner margin so pixels line up
            L = [zeros(settings.imageInnerMargin, size(imgPlot,2));
                zeros(size(L,1), settings.imageInnerMargin), L, zeros(size(L,1), settings.imageInnerMargin);
                zeros(settings.imageInnerMargin, size(imgPlot,2))];
            
            %Filter a bit
            L = movmean(movmean(L~=0,settings.pixelAverageWidth),settings.pixelAverageWidth,2);
            
            %Threshold into binary image
            L = L-(1-settings.pixelAverageThreshold)>0;
            
            timesChanged = sum(diff(L,1,1)>0,1);
            
            ximg = round(linspace(1, size(imgPlot,2), settings.nXpoints))';
            
            
            %% convolution with a triangular pulse 
            % to find middle of largest section of line as start point
            % --helps to avoid noise and legends
            preInt = round(size(L,2)/2);
            h = 1:preInt;
            h = [h, flip(h)]/max(h);

            goodLoc = 1*(timesChanged'==nTraces);


            % Seems like sometimes this is a row, sometimes a column...
            if (size(goodLoc,1) == 1) && (size(h,1) == 1)
                filtFind = conv(h,goodLoc);
                filtFind = filtFind(preInt:end-preInt).*goodLoc;
            else
                filtFind = conv(h,goodLoc');
                filtFind = filtFind(preInt:end-preInt).*goodLoc';
            end

            % start from wherever we're most confident we've found the
            % line
            startLoc = find(filtFind(ximg) == max(filtFind(ximg)),1);

            
            xloc = ximg(startLoc);
            
            % Find verticle centroids of the nTraces lines at xloc
            
            if nTracesNF > nFound
                startPtsY = obj.findClusterCentroids(L(:,xloc)~=0, 200);
            else
                startPtsY = obj.findClusterCentroids(L(:,xloc)~=0);
            end
%             preIntY = round(size(L,1)/50);
%             hY = 1:preIntY;
%             hY = [hY, flip(hY)]/max(hY);
%             
%             traceFind = conv(hY,L(:,xloc)~=0);
%             traceFind = traceFind(preIntY:end-preIntY);
%             
%             startPtsY = find(traceFind > max(traceFind)*.5 & ...
%                 [0; diff(traceFind)] >= 0 & ...
%                 [diff(traceFind); 0] <= 0);
%             
%             startPtsY([inf; diff(startPtsY)] < 2) = []; % remove any back-to-back points

            if(0)
                %% Not currently used, but finds all the points -- maybe can help with gaps due to labels/legend
                yVals = nan*ones(max(timesChanged) ,settings.nXpoints);
                for i = 1:length(ximg)
                    lineYlocs = obj.findClusterCentroids(L(:,ximg(i))~=0, 200);
                    yVals(1:length(lineYlocs), i) = lineYlocs;
                end
            end
            
            if(0)
                %% Debug testing
%                 imshow(L);
%                 hold on;
%                 plot([xloc,xloc], [1 size(imgPlot,2)],':r');
%                 plot(xloc*ones(size(startPtsY)), startPtsY, 'ob');
%                 
%                 if exist('yVals', 'var')
%                     for i=1:size(yVals,2)
%                         plot(ximg(i), yVals(:,i), 'ob');
%                     end
%                 end
%                 
%                 if exist('trace', 'var')
%                     for i = 1:size(trace,2)
%                         plot(ximg, trace(:,i), ':g', 'LineWidth',3);
%                     end
%                 end
            end
            
            trace = nan*ones(50,nTraces);
            for i = 1:length(startPtsY)
                yimg = nan*ones(1, settings.nXpoints);
                yimg(startLoc) = startPtsY(i);
                 %% going right
                for j = startLoc+1:settings.nXpoints-1
                    xloc = ximg(j);
                    Ypts = find(L(:,xloc) == 1);
                    
                    if j==startLoc+1
                        predY = yimg(j-1);
                        Yopts = Ypts;
                        Yopts([inf; diff(Yopts)] < 2) = [];
                        if length(Yopts) == length(startPtsY)
                            predY = Yopts(i);
                        end
                    else
                        predY = 2*yimg(j-1) - yimg(j-2);
                    end
                    yjump = abs(Ypts-predY);
                    [yjump, closeInd] = min(yjump);
                    onSameLine = obj.findConnectedPoints(Ypts, closeInd);
                    Ypts = Ypts(onSameLine==1);
                    
                    if ~isempty(Ypts) && yjump < size(imgPlot,1)*settings.maxYRelPtJump
                        yimg(j) = Ypts(ceil(length(Ypts)/2));
                    end

                end
                
                 %% going left
                for j = flip(1:startLoc-1)
                    xloc = ximg(j);
                    Ypts = find(L(:,xloc) == 1);
                    
                    if j==startLoc-1
                        predY = yimg(j+1);
                        Yopts = Ypts;
                        Yopts([inf; diff(Yopts)] < 2) = [];
                        if length(Yopts) == length(startPtsY)
                            predY = Yopts(i);
                        end
                    else
                        predY = 2*yimg(j+1) - yimg(j+2);
                    end
                    yjump = abs(Ypts-predY);
                    [yjump, closeInd] = min(yjump);
                    onSameLine = obj.findConnectedPoints(Ypts, closeInd);
                    Ypts = Ypts(onSameLine==1);
                    
                    if ~isempty(Ypts) && yjump < size(imgPlot,1)*settings.maxYRelPtJump
                        yimg(j) = Ypts(ceil(length(Ypts)/2));
                    end

                end
                
            trace(:,i) = yimg;
            end
            if ~isempty(trace)
                obj.loadAutoColorDigitizedData(ximg, trace);
            end
               
        end

    end
    methods (Static)
        
       function onSameLine = findConnectedPoints(ptIndex, targetIndex)
            %Helper function to find clusters of points that are connected
            % ptIndex is a vector of indices of points which are on the
            % lines
            % targetIndex is a reference for which cluster we want to find
            onSameLine = zeros(size(ptIndex));
            onSameLine(targetIndex) =1;
            for jj = targetIndex+1:length(ptIndex)
               onSameLine(jj) = sum(abs(diff(ptIndex((targetIndex:jj))))) == (jj-targetIndex);
            end
            for jj = 1:targetIndex-1
               onSameLine(jj) = sum(abs(diff(ptIndex((jj:targetIndex))))) == (targetIndex-jj);
            end
       end 
       
       function locs = findClusterCentroids(vec, filtRatio)
           % Helper function to locate centroids of distinct line crossings in a
           % vector (image slice). 
           if nargin == 1
               filtRatio = 50;
           end
            preIntY = round(length(vec)/filtRatio);
            hY = 1:preIntY;
            hY = [hY, flip(hY)]/max(hY);
            
            traceFind = conv(hY,vec);
            traceFind = traceFind(preIntY:end-preIntY);
            
            locs = find(traceFind > max(traceFind)*.5 & ...
                [0; diff(traceFind)] >= 0 & ...
                [diff(traceFind); 0] <= 0);
            
            locs([inf; diff(locs)] < 2) = []; % remove any back-to-back points 
       end
        
    end
end

