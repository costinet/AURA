function loadModel(obj, fn, swseq, force)
            % loadModel(obj, fn, swseq, force) (PLECS Implementation)
            % fn = path to plecs model
            % swseq = binary vector of switching states to be parsed
            % force = binary scalar.  If force == 0, the action won't
            % repeat on subsequent calls unless the file has been updated
            % or a new swseq is included.
            % REquired behavior:
            % populates fields As, Bs, Cs, Ds, Is, K
            %     switchLabels, stateLabels, outputLabels, inputLabels
            try
                plecs('set',fn,'EnableStateSpaceSplitting', 'off');
            catch
                warning('Unable to turn of State Space Splitting in PLECS.  Resulting matrices may not be complete')
            end
            
            if nargin == 3
                force = 0;
            end
            
            file = dir([fn(1:strfind(fn,'/')-1) '.slx']);
            if isempty(file)
                try
                    [ST,~] = dbstack;
                    fullpath = which([ST(end).name '.m']);
                    slashLocs = strfind(fullpath,'\');
                    file = dir([fullpath(1:slashLocs(end)) fn(1:strfind(fn,'/')-1) '.slx']);
                catch
                    file = [];
                end
            end

            if isempty(swseq)
                ssOrder = plecs('get', fn, 'StateSpaceOrder');
                swseq = zeros(1,length(ssOrder.Switches));
            end
            
            
            if ~force
                if ~isempty(obj.sourcefn)
                    if strcmp(obj.sourcefn, fn)
                        if ~isempty(file) && strcmp(obj.sourcefdate, file.date)
                            if nargin > 2
                                if all(ismember(swseq, obj.topology.swseq, 'rows'))
                                    return
                                else
                                    %% only new swseqs, just add those
                                    newSwSeq = setdiff(unique(swseq,'rows'), obj.topology.swseq, 'rows');
                                    for i = 1:size(newSwSeq,1)
                                        
                                        plecs('set', fn, 'SwitchVector', newSwSeq(i,:));
                                        names = plecs('get', fn, 'Topology');
                                        obj.topology.As(:,:,end+1) = names.A;
                                        obj.topology.Bs(:,:,end+1) = names.B;
                                        obj.topology.Cs(:,:,end+1) = names.C;
                                        obj.topology.Ds(:,:,end+1) = names.D;
                                        obj.topology.Is(:,:,end+1) = names.I;
                                        obj.topology.swseq = [obj.topology.swseq; newSwSeq(i,:)];
                                    end
                                    return
                                    
                                end
                            else
                                return
                            end
                        end
                    end
                end
            end  
            
            open_system(fn(1:strfind(fn,'/')-1),'loadonly');
            ssOrder = plecs('get', fn, 'StateSpaceOrder');

%             obj.topology.switchLabels = [ ...
%                 cellfun(@(x) x(strfind(x,'FET')+3:end), ssOrder.Switches, 'un', 0)' ...
%                 cellfun(@(x) x(strfind(x,'D')+1:end), ssOrder.Switches, 'un', 0)'
%                 ];
            obj.topology.switchLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Switches, 'un', 0);
            obj.topology.stateLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.States, 'un', 0);
            obj.topology.outputLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Outputs, 'un', 0);
            obj.topology.inputLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Inputs, 'un', 0);
            
           
            if nargin > 2
                if(force)
                    obj.topology.swseq = unique(swseq,'rows');
                    obj.topology.As = [];
                    obj.topology.Bs = [];
                    obj.topology.Cs = [];
                    obj.topology.Ds = [];
                    obj.topology.Is = [];
                    for i = 1:size(obj.topology.swseq,1)
                        plecs('set', fn, 'SwitchVector', obj.topology.swseq(i,:));
                        names = plecs('get', fn, 'Topology');
                        obj.topology.As(:,:,i) = names.A;
                        obj.topology.Bs(:,:,i) = names.B;
                        obj.topology.Cs(:,:,i) = names.C;
                        obj.topology.Ds(:,:,i) = names.D;
                        obj.topology.Is(:,:,i) = names.I;
                    end
                else
                    if isempty(obj.topology.swseq)
                        obj.topology.swseq = unique(swseq,'rows');
                        newSeq = obj.topology.swseq;
                    else 
                        newSeq = setdiff(swseq,obj.topology.swseq,'rows');
                        obj.topology.swseq = [obj.topology.swseq; newSeq];
                    end
                    if ~isempty(newSeq)
                        for i = 1:size(obj.topology.swseq,1)
                            plecs('set', fn, 'SwitchVector', obj.topology.swseq(i,:));
                            names = plecs('get', fn, 'Topology');
                            obj.topology.As(:,:,i) = names.A;
                            obj.topology.Bs(:,:,i) = names.B;
                            obj.topology.Cs(:,:,i) = names.C;
                            obj.topology.Ds(:,:,i) = names.D;
                            obj.topology.Is(:,:,i) = names.I;
                        end
                    end
                end

            end
            
            obj.topology.K = eye(length(obj.topology.stateLabels));
            
            for i = 1:length(obj.topology.stateLabels)
                cloc = strfind(obj.topology.stateLabels{i}, ':');
                if isempty(cloc)
                    element = plecs('get',[fn, '/', obj.topology.stateLabels{i}]);
                else
                    element = plecs('get',[fn, '/', obj.topology.stateLabels{i}(1:cloc-1)]);
                end
                if strcmp(element.Type, 'Capacitor')
                    obj.topology.K(i,i) = evalin('base',element.C);
                elseif strcmp(element.Type, 'Inductor')
                    obj.topology.K(i,i)  = evalin('base',element.L);
                elseif strcmp(element.Type, 'Transformer')
                    obj.topology.K(i,i)  = evalin('base',element.Lm);
                end
            end
            
            obj.sourcefn = fn;
            if isempty(file)
                obj.sourcefdate = '0';
            else
                obj.sourcefdate = file.date;
            end

        end