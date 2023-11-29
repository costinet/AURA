classdef netConnect < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        serverURL
        apiKey
        post_max_size 
        upload_max_filesize
    end
    
    methods
        function obj = netConnect(apiKey, serverURL)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin <2 
                if exist('webSettings.mat', 'file')
                    load('webSettings.mat', 'apiKey', 'serverURL');
                end
                if ~exist('apiKey', 'var') || isempty(apiKey)       
                    warning('no Database or apiKey specified, please provide them now');
                    serverURL = input('Enter URL of database server\n', 's');
                    apiKey = input('Enter apiKey \n', 's');
                end
            end
            
            obj.serverURL = serverURL;
            obj.apiKey = apiKey;
            try
                obj.handshakeServer();
            catch e
                if strcmp(e.identifier, 'MATLAB:webservices:Timeout')
                    error('Unable to initialize net connection.  No response from server');
                else
                    throw(e);
                end
            end
        end
        
        function handshakeServer(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            options = weboptions('Timeout', 10);
            data = webread([obj.serverURL '/upload?apiKey=' obj.apiKey '&action=handshake'], options);
            if strcmp(data(1:5), 'Error')
                error(data);
            end
            data = strrep(data, 'M', 'e6');
            data = str2double(strsplit(data,', '));

            valid = (length(data) == 2) && all(~isnan(data));
            assert(valid, 'Handshaking with server failed');

            obj.post_max_size = data(1);
            obj.upload_max_filesize = data(2);
        end
        
        function result = postData(obj, objDB, startIndex)
            if nargin == 2
                startIndex = 1;
            end
            
            [~, name] = system('hostname');
            
            result = '';
            
            %encodedDB = jsonencode(objDB(startIndex:end));
            newEntries = objDB(startIndex:end);
            newEntries = newEntries([newEntries.upDated] == 1);
            if isempty(newEntries)
                result = '';
                return
            end
            encodedDB = jsonencode(newEntries);

            jsonData = whos('encodedDB');

            
            if jsonData.bytes < obj.post_max_size % Size OK
                PostName1 = 'data';
                PostValue1 = encodedDB;
                result = webwrite([obj.serverURL '/upload?apiKey=' obj.apiKey '&action=upload'],PostName1,PostValue1, 'name', name);
            else
                %DB size is too big for a single POST, split up and call
                %recursively
                for nI = startIndex:length(objDB)
                    encodedDB = jsonencode(objDB(startIndex:nI));
                    jsonData = whos('encodedDB');
                    if jsonData.bytes > obj.post_max_size
                        nI=nI-1;
                        break;
                    end
                    
                end
                if nI<startIndex
                    error('transistor too big!');
                else
                    encodedDB = jsonencode(objDB(startIndex:nI));
                    PostName1 = 'data';
                    PostValue1 = encodedDB;
                    response = webwrite([obj.serverURL '/upload?apiKey=' obj.apiKey '&action=upload'],PostName1,PostValue1, 'name', name);
                    result = [result response];
                    if nI < length(objDB)
                        %Recursive Call
                        response = obj.postData(objDB, nI+1);
                        result = [result response];
                    end
                end
            end
        end
        
        function [graphs, params] = pullData(obj)
            newData = webread([obj.serverURL '/upload?apiKey=' obj.apiKey '&action=download']);
            
%             fid = fopen('C:\Users\dcostine\Desktop\new 10.html', 'w');
%             fwrite(fid, newData);
%             fclose(fid);
%             
            graphs = [];
            params = [];
            graphInd = [regexpi(newData, '<graphData>') regexpi(newData, '<\\graphData>')];
            if ~isempty(graphInd)
                newGraphs = newData(graphInd(1)+11:graphInd(end)-1);
                graphs = jsondecode(newGraphs);
            end
            
            paramInd = [regexpi(newData, '<paramData>') regexpi(newData, '<\\paramData>')];
            if ~isempty(paramInd)
                newParams = newData(paramInd(1)+11:paramInd(end)-1);
                params = jsondecode(newParams);
            end
        end
    end
end

