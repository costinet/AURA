classdef FET < handle
    %FET contains all of the data of a Field Effect Transistor (FET)
    %   Now looking at input must be of the form 'property' value then
    %   numerical value - similar to the constructors for figures
    
    % Will have to determine what variable names will be and how we
    % can change them based on availible data
    
%     %%%%%%   %      %  %%%%%%%    %%%%%%
%    %      %  %      %  %      %  %      %
%    %      %  %      %  %      %  %      %
%    %%%%%%%%  %      %  %%%%%%%   %%%%%%%%
%    %      %  %      %  %%        %      %
%    %      %  %      %  % %       %      %
%    %      %  %      %  %  %      %      %
%    %      %  %      %  %   %     %      %
%    %      %   %    %   %    %    %      %
%    %      %    %%%%    %     %   %      %
    
    properties

        R_on
        R_diode
        V_f
        R_off
        
    end

    methods

        function initialize(~,nargin)
            % This function creates variables for all of the FET
            % properties
            
            if mod(nargin,2)~=0
                fprintf('input length is not correct \n')
            end
            
            
            
            for i = 1:1:length(nargin)
            
                
                if mod(i,2)
                    % This should be the string identifier for the
                    % next 
                    
                    if ~ischar(nargin(i))
                        fprintf('Not correct class type for odd entries')
                        
                    end
                    
                else
                    % This should be the value of the previous string identifier
                    
                    if ~isnumeric(nargin(i))
                        fprintf('Not correct class type for even entries')
                    end
                end
            end
        end

    end
    % That's all Folks
end
