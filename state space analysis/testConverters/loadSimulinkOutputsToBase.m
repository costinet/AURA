function loadSimulinkOutputsToBase(modelfile,PLECsBlock)
%loadSimukinkOutputsToBase load all variables in simout into the base
%workspace
%
%   This is a helper function for easy testing of Simulink-embedded PLECS
%   files such as the example topologies provided here.  The function
%   simulates the simulink model for a minimum time and then transfers all
%   variables loaded into simout into the base workspace so they are
%   available for PLECS to read
%
%   Note that this WILL overwrite variables in the base workspace of the
%   same name, so use with care.

    open_system(modelfile,'loadonly');
    circuitPath = [modelfile '/' PLECsBlock];
    set_param(circuitPath,'Commented','on');
    simout = sim(modelfile,eps);
    
    for i = 1:length(simout.properties)
        assignin('base',simout.properties{i},eval(['simout.' simout.properties{i}]));
    end
    
    set_param(circuitPath,'Commented','off');
end