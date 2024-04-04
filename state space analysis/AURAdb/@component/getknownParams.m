function getknownParams()
        % getknownParams() is used to convert the csv file into text strings that 
        % can be copied into the constant variables
        [~,~,raw] = xlsread('KnownParams.csv');
        paramNames = '{';
        knownParams = '{';
        defaultUnits = '{';
        defaultMultipliers = '{';
        dictKeys = '{';
        dictVals = '{';
        for i = 1:size(raw,2)
            paramNames = [paramNames, '''', raw{1,i},''','];
            knownParams = [knownParams, '''', raw{2,i},''','];
            if isletter(raw{3,i}(1))
                defaultUnits = [defaultUnits, '''', raw{3,i},''','];
            elseif strcmp(raw{3,i}(1), '[')
                loc = strfind(raw{3,i},']');
                defaultUnits = [defaultUnits, '''', eval(raw{3,i}(1:loc)), raw{3,i}(loc+1:end),''','];
            end
            
            if isletter(raw{4,i}(1))
                defaultMultipliers = [defaultMultipliers, '''', raw{4,i},''','];
            elseif isnumeric(raw{4,i}(1))
                if raw{4,i}(1) == 1
                    defaultMultipliers = [defaultMultipliers, '''',''','];
                else
                    error('invalid multiplier');
                end
            elseif strcmp(raw{4,i}(1), '[')
                loc = strfind(raw{4,i},']');
                defaultMultipliers = [defaultMultipliers, '''', eval(raw{4,i}(1:loc)), raw{4,i}(loc+1:end),''','];
            end
%             defaultMultipliers = [defaultMultipliers, '''', num2str(raw{4,i}),''','];
            for j = 5:size(raw,1)
                if ~isnan(raw{j,i})
                    dictKeys = [dictKeys, '''', raw{j,i},''','];
                    dictVals = [dictVals, '''', raw{1,i},''','];
                end
            end
        end
        paramNames = [paramNames(1:end-1), '}'];
        knownParams = [knownParams(1:end-1), '}'];
        defaultUnits = [defaultUnits(1:end-1), '}'];
        defaultMultipliers = [defaultMultipliers(1:end-1), '}'];
        paramDict = ['containers.Map(' dictKeys(1:end-1) '},' dictVals(1:end-1) '})'];

        paramDef = ['knownParams = ' knownParams newline ....
            'defaultUnits = ' defaultUnits newline ...
            'defaultMultipliers = ' defaultMultipliers newline...
            'paramDict = ' paramDict newline ...
            'paramNames = ' paramNames];

        disp(paramDef)
        
%         SQLtable = 

    end
