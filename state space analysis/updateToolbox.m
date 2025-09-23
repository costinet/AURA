%% updateToolbox
% This script checks the current version of the Switched Mode Power Supply
% Toolbox against the latest github release.  If the currently installed
% version predates the latest release, it will download and install the
% newer version

gitRelease = webread('https://api.github.com/repos/costinet/aura/releases/latest');
latestVersion = gitRelease.tag_name;
% installedVersion = ver("Switched Mode Power Supply Toolbox").Version;
% isMATLABReleaseOlderThan('R2023b')
toolboxes = matlab.addons.toolbox.installedToolboxes;
if isempty(toolboxes)
    thisToolbox = {};
else
    thisToolbox = toolboxes(strcmp({toolboxes.Name},"Switched Mode Power Supply Toolbox"));
end
if isempty(thisToolbox)
    releaseURL = gitRelease.assets.browser_download_url;
    [~, defaultName, ext] = fileparts(gitRelease.assets.name);
    outfn = fullfile(tempdir, [defaultName, ext]);
    
    selection = questdlg(['Found no existing installation of the Switched Mode Power Supply Toolbox. ' ...
        'The latest release from github is ' latestVersion '. Would you like to install?'], ...
        "Confirm Install");
    if strcmp(selection, 'Yes')
        websave(outfn, releaseURL);
        installedToolbox = matlab.addons.toolbox.installToolbox(outfn);
        delete(outfn);

        % Update Component Libraries
        AURAdb(1).updateLibraries;
    end
    return
else
    installedVersion = thisToolbox.Version;
end

% verDiff = cellfun(@str2num,split(latestVersion(2:end),'.')) - cellfun(@str2num,split(installedVersion,'.'));
% if ~isempty(find(verDiff > 0,1,'first') )
%     if isempty(find(verDiff < 0,1,'first')) || ...
%         find(verDiff > 0,1,'first') < find(verDiff < 0,1,'first')
% if verLessThan(installedVersion,latestVersion(2:end)) % ver() not working
% in latest releases
if compareReleases(latestVersion, installedVersion)
        releaseURL = gitRelease.assets.browser_download_url;
        [~, defaultName, ext] = fileparts(gitRelease.assets.name);
        outfn = fullfile(tempdir, [defaultName, ext]);

        websave(outfn, releaseURL);

         selection = questdlg(['Found version ' installedVersion ' of the Switched Mode Power Supply Toolbox installed. ' ...
        'The latest release from github is ' latestVersion '. Would you like to install the latest version?'], ...
        "Confirm Install");
        if strcmp(selection, 'Yes')
            matlab.addons.toolbox.uninstallToolbox(thisToolbox);
            installedToolbox = matlab.addons.toolbox.installToolbox(outfn);
            delete(outfn);
        end
    % end
else
    disp('Toolbox is up to date')
end

%% update component libraries
AURAdb(1).updateLibraries

function result = compareReleases(latestVersion, installedVersion)
    toolboxParts = getParts(installedVersion);
    verParts = getParts(latestVersion);
    if toolboxParts(1) ~= verParts(1)     % major version
        result = toolboxParts(1) < verParts(1);
    elseif toolboxParts(2) ~= verParts(2) % minor version
        result = toolboxParts(2) < verParts(2);
    else                                  % revision version
        result = toolboxParts(3) < verParts(3);
    end
end

function parts = getParts(V)
    % parts = sscanf(V, '%d.%d.%d')';
    tokens = str2double(regexp(V, 'v?(\d+)\.(\d+)\.(\d+)', 'tokens', 'once'));
    if length(parts) < 3
        parts(3) = 0; % zero-fills to 3 elements
    end
end


