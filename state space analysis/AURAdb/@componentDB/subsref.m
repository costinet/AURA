function varargout = subsref(obj,s)
    [varargout{1:nargout}] = builtin('subsref',obj,s);
    
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
end

