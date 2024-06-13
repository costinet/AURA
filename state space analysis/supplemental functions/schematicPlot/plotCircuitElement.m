function plotCircuitElement(f,type, loc, sz, rot, baseLineWidth)
    arguments
        f
        type {mustBeTextScalar,mustBeMember(type,{'M','C','L','R','V','I','D','Vm'})}
        loc (1,2) {mustBeNumeric}
        sz (1,1) {mustBeNumeric, mustBeScalarOrEmpty}
        rot (1,1) {mustBeNumeric, mustBeScalarOrEmpty}
        baseLineWidth (1,1) {mustBeNumeric, mustBeScalarOrEmpty}   
    end
   
    if isempty(f)
        f = figure;
    else
        assert(isa(f,'matlab.ui.Figure'), 'Input f must be a handle to a figure')
    end

    [skinnyLines,thickLines,fills] = circuitElementImage(type);
   
    skinnyLines = scaleAndRot(skinnyLines, sz, rot, loc);
    thickLines = scaleAndRot(thickLines, sz, rot, loc) ;
    fills = scaleAndRot(fills, sz, rot, loc) ;
    % 
    figure(f); hold on;
    plot(skinnyLines(:,1),skinnyLines(:,2),'k','LineWidth',baseLineWidth);
    if ~isempty(thickLines)
        plot(thickLines(:,1),thickLines(:,2),'k','LineWidth',2*baseLineWidth);
    end
    if ~isempty(fills)
        fill(fills(:,1),fills(:,2),'k','LineWidth',baseLineWidth);
    end


end

function newPts = scaleAndRot(points, sz, rot, loc) 
    if isempty(points)
        newPts = [];
        return
    end
    points = points*sz/10;
    th = atan2(points(:,2),points(:,1)) + rot/180*pi;
    points = sqrt(sum(points.^2,2)).*[cos(th) sin(th)];
    points = points + repmat([loc(1) loc(2)],size(points,1),1);
    newPts = points;
end