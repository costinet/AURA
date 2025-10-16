function smps = loadPLECSTestModel(modelfile, PLECsModel)
%smps = loadPLECSTestModel(modelfile, PLECsModel) loads a formatted
%simulink file containing a test circuit for the switched mode power supply
%toolbox
%   This function calls loadSimulinkOutputsToBase and implements error
%   catching for improperly-formatted test setups and for users without
%   PLECS installed.  It is primarily used in the example files to show
%   PLECS operation, but allow users to execute the code even if PLECS is
%   not installed.
%   
%   see also SMPSim, loadSimulinkOutputsToBase, plecs

    try
        loadSimulinkOutputsToBase(modelfile,PLECsModel);
        circuitPath = [modelfile '/' PLECsModel];
        smps = SMPSim();
        swvec = evalin('base','swvec');
        us = evalin('base','us');
        ts = evalin('base','ts');
        smps.initialize(circuitPath, swvec, us, ts);   
    catch e
        if strcmp(e.identifier, 'LOADCIRCUIT:plecsNotInstalled')
            warning(e.message)
            smps = SMPSim();
        elseif strcmp(e.identifier, 'MATLAB:UndefinedFunction')
            warning('Simulink model did not define one of the required variables: swvec, us, ts.  The returned simulation object will not be runnable');
        else 
            rethrow(e)
        end
    end

end