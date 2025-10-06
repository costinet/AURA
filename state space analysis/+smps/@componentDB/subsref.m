function varargout = subsref(obj, s)
    % Overload subsref for the componentDB class.
    
%     switch s(1).type
%         case '()'
%             % --- Handle Parenthesis Indexing: db(1:3) or db(1).name ---
%             try
%                 indexed_components = obj.components(s(1).subs{:});
%             catch e
%                 if isempty(obj.components)
%                     ME = MException('componentDB:empty', 'The componentDB database is empty. Cannot index into it.');
%                     e = e.addCause(ME);
%                 end
%                 rethrow(e);
%             end
% 
%             if length(s) == 1
%                 % Simple indexing, e.g., db(5). Result is the component object(s).
%                 varargout = {indexed_components};
%             else
%                 % This is chained indexing, e.g., db(1:3).name or [v,t]=db(1:3).myMethod()
% 
%                 num_components = numel(indexed_components);
%                 % The number of outputs the user wants from the entire expression.
%                 num_outputs_requested = nargout;
% 
%                 % Pre-allocate the final output cell arrays.
%                 % Each varargout{j} will hold all the j-th outputs from all components.
%                 varargout = cell(1, num_outputs_requested);
%                 for j = 1:num_outputs_requested
%                     varargout{j} = cell(1, num_components);
%                 end
% 
%                 % Loop through each component that was indexed.
%                 for i = 1:num_components
%                     % Create a temporary cell array to capture the outputs from this
%                     % single component's subsref call.
%                     temp_outputs = cell(1, num_outputs_requested);
% 
%                     try
%                         % This is the key step. We ASK the inner subsref for the number
%                         % of outputs the user wants.
%                         [temp_outputs{:}] = subsref(indexed_components(i), s(2:end));
%                     catch e
%                         % If the above line fails, it's most likely because we asked
%                         % for too many outputs (e.g., asking .name for two outputs).
%                         if strcmp(e.identifier, 'MATLAB:TooManyOutputs')
%                             % This is not a bug to hide, it's a user error.
%                             % Throw a new, more informative error.
%                             ME = MException('componentDB:TooManyOutputs', ...
%                                 'The chained expression cannot produce %d outputs. Check property names and method signatures.', ...
%                                 num_outputs_requested);
%                             throw(ME);
%                         else
%                             % It was some other legitimate error.
%                             rethrow(e);
%                         end
%                     end
% 
%                     % If the call succeeded, distribute the results into the final
%                     % output cell arrays.
%                     for j = 1:num_outputs_requested
%                         varargout{j}{i} = temp_outputs{j};
%                     end
%                 end
%             end
% 
%         case '.'
%             % --- Handle Dot Indexing: db.description OR db.value ---
%             try
%                 [varargout{1:nargout}] = builtin('subsref', obj, s);
%             catch e
%                 if strcmp(e.identifier, 'MATLAB:noSuchMethodOrField')
%                     if isempty(obj.components)
%                         rethrow(e);
%                     end
%                     [raw_output{1:nargout}] = arrayfun(@(comp) subsref(comp, s), ...
%                                                       obj.components, 'UniformOutput', false);
% 
%                     for i = 1:nargout
%                        if all(cellfun(@(x) isnumeric(x) && isscalar(x), raw_output{i}))
%                            varargout{i} = [raw_output{i}{:}];
%                        else
%                            varargout{i} = raw_output{i};
%                        end
%                     end
%                 else
%                     rethrow(e);
%                 end
%             end
% 
%         case '{}'
%             error('componentDB:badsubs', 'Brace indexing {} is not supported for componentDB objects.');
%     end
% end

    try 
        

        % One case we don't want the standard behavior is if we index the
        % first element
        if strcmp(s(1).type, '()') && s(1).subs{1} == 1
            if isscalar(s)
                varargout = {obj.components(1)};
            else
                [varargout{1:nargout}] = builtin('subsref',obj.components(1),s(2:end));
            end
        elseif strcmp(s(1).type, '()') && strcmp(s(1).subs{1},':')
            if isscalar(s)
                varargout = {obj.components(:)};
            else
                [varargout{1:nargout}] = builtin('subsref',obj.components(:),s(2:end));
            end
        else
            %prioritize standard behavior
            [varargout{1:nargout}] = builtin('subsref',obj,s);
        end
        return
    catch e
        try
             %numerical indexing is passed on to components list
            if strcmp(s(1).type, '()')
                if isscalar(s)
                    if isnumeric(s(1).subs{:}) && numel(obj.components) >=  max(s(1).subs{:})
                        varargout = {obj.components(s(1).subs{:})};
                    elseif ischar(s(1).subs{:}) && strcmp(s(1).subs{:},':')
                        varargout = {obj.components(:)};
                    else
                        try 
                            varargout = {obj(s(1).subs{:})};
                        catch e3
                            if isempty(obj)
                                ME = MException('componentDB:empty', ...
                                    ['AURAdb databse ' class(obj) ' is empty.  Run updateToolbox from the command line to get the latest databases from the repository']);
                                e3 = e3.addCause(ME);
                            end
                            rethrow(e3)
                        end
                    end
                else 
                    comps = obj.components(s(1).subs{:});
                    varargout = {};
                    if length(comps) == 1
                        component = comps(1);
                        varargout = {subsref(component,s(2:end))};
                        % actNargout = numel(result);
                        % [varargout{1:actNargout}] = result(:);
                    else
                        for i = 1:length(comps)
            %                 [varargout{i}] = builtin('subsref',comps(i),s(2:end));
                            component = comps(i);
                            res = subsref(component,s(2:end));
                            if ischar(res)
                                res = cellstr(res);
                            end
                            if numel(res) > 1
                                [compOut{1:nargout}] = subsref(component,s(2:end));
                                [varargout{i,1:nargout}] = compOut{1:nargout};
                            else
                                [varargout{i}] = subsref(component,s(2:end));
                            end
                        end
                        for i = 1:nargout
                            varargout{i} = {varargout{:,i}};
                        end
                    end
                end        
            else
                try
                    %try applying it to the DB
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                catch e
                    %If that fails, try applying it to the components in the DB
                    try
                        [varargout{1:nargout}] = arrayfun(@(x)builtin('subsref',x,s),obj.components,'UniformOutput',false); 
                    catch e2
                        % If that fails, throw the original error from trying to
                        % apply it to the component
                        throw(e);
                    end
                end
            end



        catch e2
            e.addCause(e2)
            rethrow(e)
        end

    end


%     [varargout{1:nargout}] = builtin('subsref',obj,s);
   
end
    
    %% MATLAB TEMPLATE
%    switch s(1).type
%       case '.'
%          if length(s) == 1
%             % Implement obj.PropertyName
%             ...
%          elseif length(s) == 2 && strcmp(s(2).type,'()')
%             % Implement obj.PropertyName(indices)
%             ...
%          else
%             [varargout{1:nargout}] = builtin('subsref',obj,s);
%          end
%       case '()'
%          if length(s) == 1
%             % Implement obj(indices)
%             ...
%          elseif length(s) == 2 && strcmp(s(2).type,'.')
%             % Implement obj(ind).PropertyName
%             ...
%          elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
%             % Implement obj(indices).PropertyName(indices)
%             ...
%          else
%             % Use built-in for any other expression
%             [varargout{1:nargout}] = builtin('subsref',obj,s);
%          end
%       case '{}'
%          if length(s) == 1
%             % Implement obj{indices}
%             ...
%          elseif length(s) == 2 && strcmp(s(2).type,'.')
%             % Implement obj{indices}.PropertyName
%             ...
%          else
%             % Use built-in for any other expression
%             [varargout{1:nargout}] = builtin('subsref',obj,s);
%          end
%       otherwise
%          error('Not a valid indexing expression')
%    end

%% OLD ATTEMPT

%  function varargout = subsref(obj, param)
% % %     try
% % %         [varargout{1:nargout}] = builtin('subsref',obj,param);
% % %     catch bultInError
% % %         try 
% % %            nargout = length(obj.components);
% % %            [varargout{1:nargout}] = builtin('subsref',obj.components,param);
% % %         catch
% % %             throw(bultInError)
% % %         end
% % %     end
% 
% 
%     % Overload dot-indexing
%     if length(param) == 1
%         if strcmp(param.type, '()')
%             if all(floor([param.subs{:}]) == [param.subs{:}]) % all integers
%                 varargout{1} = obj.components([param.subs{:}]);
%             end
%         else
%             [varargout{1:nargout}] = builtin('subsref',obj,param);
%         end
%     elseif length(param) == 2
%         if strcmp(param(1).type, '()')
%             if all(floor([param(1).subs{:}]) == [param(1).subs{:}]) % all integers
%                 [varargout{1:nargout}] = builtin('subsref',obj.components([param(1).subs{:}]),param(2));
%             end
%         end
%     else
%        [varargout{1:nargout}] = builtin('subsref',obj,param); 
%     end
% end

