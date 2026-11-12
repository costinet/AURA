function normErr = compareResultToLTSpice(sim, figNo, options)
%COMPARERESULTTOLTSPICE compares toolbox steady-state prediction to LTspice
%simulation
%   normErr = COMPARERESULTTOLTSPICE(sim, figNo, leaveNetlist)
arguments (Input)
    sim SMPSim
    figNo {isinteger}
    options.leaveNetlist (1,1) {islogical} = 0
    options.fileIOTime (1,1) {isinteger} = 2
end

% arguments (Output)
%     normErr
% end

assert(isa(sim.parser, 'NetlistCircuitParser'), 'Validation with LTSpice requires the use of NetlistCircuitParser');

nl = sim.parser.createSimulateableNetlist(sim);
fn = 'temp';
if exist([fn '.net'],'file')
    i=2;
    while exist(sprintf([fn '%i.net'], i),'file')
        i=i+1;
    end
    fn = sprintf([fn '%i'], i);
else
    fn = 'temp';
end



fid = fopen([fn '.net'], 'w');
fwrite(fid, nl);
fclose(fid);



try
    assert(waitForFile(fullfile(pwd, [fn '.net']), options.fileIOTime), 'File not written in time')

    result = system([ '"'...
        sim.parser.LTSpiceExe '"' ' -b -run "'...
        fullfile(pwd, [fn '.net']) '"']);
    
    assert(result == 0, 'Failed running LTSpice')

    assert(waitForFile(fullfile(pwd, [fn '.raw']), options.fileIOTime), 'Raw file not available in time')
    data = LTspice2Matlab2( fullfile(pwd, [fn '.raw']) );

catch e
    delete(fullfile(pwd, [fn '.net']));
    if exist(fullfile(pwd, [fn '.log']), 'file')
        fid = fopen(fullfile(pwd, [fn '.log']), 'r');
        % logText =  fread(fid, '*char')';
        data = fread(fid, '*uint8')';
        logText = char(typecast(data, 'uint16'));
        fclose(fid);
        delete(fullfile(pwd, [fn '.log']));

        log = MException('LTSpice:LoggedError', ['LTspice log information: \n' logText]);
        e = e.addCause(log);
    end
    rethrow(e)
end



if ~options.leaveNetlist
    delete(fullfile(pwd, [fn '.*']));
end

figure(figNo);
clf(figNo);
sim.plotAllStates(figNo);
fig = figure(figNo);
hold(fig.Children, 'on');

for i = 1:numel(sim.stateNames)
    
    if startsWith(sim.stateNames(i), 'L')
        curAxes = findobj(fig,'Tag',['I_{' sim.stateNames{i} '}']);
        plot(curAxes, data.time_vect,data.variable_mat( ...
            strcmpi(data.variable_name_list,['I(',sim.stateNames{i} ')'])   ,:), ...
            ':r','LineWidth',2)
    elseif startsWith(sim.stateNames(i), 'C')
        if startsWith(sim.stateNames(i), 'C_M') || startsWith(sim.stateNames(i), 'C_D')
            curAxes = findobj(fig,'Tag',['V_{C_{' sim.stateNames{i}(3:end) '}}']);
            plot(curAxes, data.time_vect,data.variable_mat( ...
                strcmpi(data.variable_name_list,['V(v_',sim.stateNames{i} ')'])   ,:), ...
                ':r','LineWidth',2)
        else
            curAxes = findobj(fig,'Tag',['V_{' sim.stateNames{i} '}']);
            plot(curAxes, data.time_vect,data.variable_mat( ...
                strcmpi(data.variable_name_list,['V(v_',sim.stateNames{i} ')'])   ,:), ...
                ':r','LineWidth',2)
        end
    end

    simData = [curAxes.Children(1).XData; curAxes.Children(1).YData];
    simData(isnan(simData)) = 0;
    spiceData = [curAxes.Children(end).XData; curAxes.Children(end).YData];

    % normErr(i) = norm(simData(2,:) - ...
    %     interp1(spiceData(1,:), spiceData(2,:), simData(1,:), 'linear', 'extrap')...
    %     ) / norm(simData(2,:) );

    normErr(i) = trapz(simData(1,:), (simData(2,:) - ...
            interp1(spiceData(1,:), spiceData(2,:), simData(1,:), 'linear', 'extrap')).^2);
    normErr(i) = normErr(i)/trapz(simData(1,:), (simData(2,:)).^2);

end


end

function result = waitForFile(path, maxTime)
    ctr = 0;
    while ~exist(path, "file") && ctr < maxTime
        pause(1);
        ctr=ctr+1;
    end
    result = ~(ctr >= maxTime);
end




