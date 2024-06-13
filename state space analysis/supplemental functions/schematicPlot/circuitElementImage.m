function [skinnyLines,thickLines,fills] = circuitElementImage(type)
%circuitElementImage returns point vectors for plotting basic circuit
%elements
%   Possible type values are 'nmos','cap','ind','res','vsrc','isrc'
    if strcmpi(type,'M')
         skinnyLines = [   -3.1250         0
           -4.1667         0
               NaN       NaN
           -3.1250   1.8750
           -3.1250   -1.8750
               NaN       NaN
                 0   5.0000
                 0   1.8750
           -1.8750   1.8750
               NaN       NaN
           -1.8750         0
                 0         0
                 0    -5.0000
               NaN       NaN
           -1.8750    -1.8750
                 0    -1.8750
               NaN       NaN
           -1.8750   3.1250
           -1.8750    -3.1250];
        
        thickLines = [   -3.1250   1.8750
           -3.1250    -1.8750];
        
        fills = [  -1.6667         0
           -1.0417   0.3604
           -1.0417    -0.3604
           -1.6667         0];
    elseif strcmpi(type,'C')
        skinnyLines = [0    5
             0    0.6250
           NaN       NaN
             0   -5
             0   -0.6250];
    
        thickLines = [-1.5625   -0.6250
            1.5625   -0.6250
               NaN       NaN
           -1.5625    0.6250
            1.5625    0.6250];

        fills = [];

    elseif strcmpi(type,'L')
        skinnyLines = [25.5, 56
            25.5, 44];
        numPtsCurves = 10;
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [-2.5,0], [-4.5,-2], [-4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [0,-2.5], [2, -4.5], [4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [.83,0], [1.5,0.67], [1.5,1.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [0,.83], [-0.67,1.5], [-1.5,1.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [-2.5,0], [-4.5,-2], [-4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [0,-2.5], [2, -4.5], [4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [.83,0], [1.5,0.67], [1.5,1.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [0,.83], [-0.67,1.5], [-1.5,1.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [-2.5,0], [-4.5,-2], [-4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [0,-2.5], [2, -4.5], [4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [.83,0], [1.5,0.67], [1.5,1.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [0,.83], [-0.67,1.5], [-1.5,1.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [-2.5,0], [-4.5,-2], [-4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; skinnyLines(end,:) + cubicBezier([0,0], [0,-2.5], [2, -4.5], [4.5,-4.5], numPtsCurves)];
        skinnyLines = [skinnyLines; 25.5, skinnyLines(end,2)-12];
        skinnyLines = [skinnyLines(:,1)-25.5, skinnyLines(:,2)-30.5]/51*10;

        thickLines = [];
        fills = [];

    elseif strcmpi(type,'R')
         skinnyLines = [0   -5.0000
                     0   -1.1538
                0.5769   -0.8654
                0.5769   -0.7692
               -0.5769   -0.1923
                0.5769    0.3846
               -0.5769    0.9615
               -0.5769    1.0577
                     0    1.3462
                     0    5.0000];
        thickLines = [];
        fills = [];
    elseif strcmpi(type,'V')
        numPtsCurves = 10;
        th = linspace(0,2*pi+pi/(2*numPtsCurves),4*numPtsCurves)';
        r = 1.7308;
        thickLines = [r*cos(th) , r*sin(th)];

        skinnyLines = [         0    1.7308
                 0    5.0000
               NaN       NaN
                 0   -1.7308
                 0   -5.0000
               NaN       NaN
                 0    1.3038
                 0    0.2577
               NaN       NaN
           -0.5231    0.7808
            0.5231    0.7808
               NaN       NaN
           -0.4327   -0.7615
            0.4212   -0.7615];

        fills = [];
     elseif strcmpi(type,'I')
        numPtsCurves = 10;
        th = linspace(0,2*pi+pi/(2*numPtsCurves),4*numPtsCurves)';
        r = 1.7308;
        thickLines = [r*cos(th) , r*sin(th)];

        skinnyLines = [  0   -1.7308
                 0   -5.0000
               NaN       NaN
                 0    1.7308
                 0    5.0000
               NaN       NaN
                 0   -1.1538
                 0    0.1923];

        fills = [ 0    0.0962
           -0.3846    0.0962
                 0    1.2500
            0.3846    0.0962
                 0    0.0962];
   elseif strcmpi(type,'D')
         skinnyLines = [0 -5
                0 5];
        thickLines = [-1   -0.75
                1   -0.75];
        fills = [0   -0.75
                1   1
               -1    1
                 0   -0.75];
    elseif strcmpi(type,'Vm')
        skinnyLines = [0 -5
                0 5];
        thickLines = [];
        fills = [0   -0.75
                .5   1
               -.5    1
                 0   -0.75];        
    else
        error('Invalid component type');
    end
end


function [pts] = cubicBezier(st, pt1, pt2, ed, numpts)
    t = linspace(0,1,numpts)';
    pts = (1-t).^3*st + 3*(1-t).^2.*t*pt1 + 3*(1-t).*t.^2*pt2 + t.^3*ed;
end