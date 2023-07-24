% post_max_size = 8e6;
% upload_max_filesize = 2e6;
% max_file_uploads = 20;

clear all 
close all
clc
addpath(genpath(pwd))

deviceNumbers = Transistor.listDevices;
for i = 1:1;%length(deviceNumbers)
    t{i} = Transistor(deviceNumbers{i}, 1);
end

URL = 'http://powerlib.ddns.net/powerlib';

%% Check capabilities
data = webread([URL '/upload?apiKey=qR06gPrxTo&action=handshake']);
if strcmp(data(1:5), 'Error')
    error(data);
end
data = strrep(data, 'M', 'e6');
data = str2double(strsplit(data,', '));

valid = (length(data) == 2) && all(~isnan(data));
assert(valid, 'Handshaking with server failed');

post_max_size = data(1);
upload_max_filesize = data(2);

% Confirm server file size limits
encodedTransistor = jsonencode(t);
jsonData = whos('encodedTransistor');
if jsonData.bytes > min(post_max_size, upload_max_filesize)
    error('requested upload too large');
end


%% Upload it
% devices = Transistor.listDevices;
[~, name] = system('hostname');


PostName1 = 'data';
PostValue1 = encodedTransistor;

RESPONSE = webwrite([URL '/upload?apiKey=qR06gPrxTo&action=upload'],PostName1,PostValue1, 'name', name)





return
error('test');
if(0)
    query = sprintf('ALTER table transistor_table\nADD ');
    lascol = 'submittedBy';
    for i = 1:length(t{1}.tableFields)
        tableColumnAddString{1} = sprintf('%s_min FLOAT,\n',t{1}.tableFields{i});
        tableColumnAddString{2} = sprintf('%s_typ FLOAT,\n',t{1}.tableFields{i});
        tableColumnAddString{3} = sprintf('%s_max FLOAT,\n',t{1}.tableFields{i});
        query = [query tableColumnAddString{1:3}]
    end
    
end