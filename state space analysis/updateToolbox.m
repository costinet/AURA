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
thisToolbox = toolboxes(strcmp(toolboxes.Name,"Switched Mode Power Supply Toolbox"));
if isempty(thisToolbox)
    releaseURL = gitRelease.assets.browser_download_url;
    outfn = websave(gitRelease.assets.name, releaseURL);
    selection = questdlg(['Found no existing installation of the Switched Mode Power Supply Toolbox. ' ...
        'The latest release from github is ' latestVersion '. Would you like to install?'], ...
        "Confirm Install");
    if strcmp(selection, 'Yes')
        installedToolbox = matlab.addons.toolbox.installToolbox(outfn);
    end
    return
else
    installedVersion = thisToolbox.Version;
end

verDiff = cellfun(@str2num,split(latestVersion(2:end),'.')) - cellfun(@str2num,split(installedVersion,'.'));
if ~isempty(find(verDiff > 0,1,'first') )
    if isempty(find(verDiff < 0,1,'first')) || ...
        find(verDiff > 0,1,'first') < find(verDiff < 0,1,'first')
        releaseURL = gitRelease.assets.browser_download_url;
        outfn = websave(gitRelease.assets.name, releaseURL);

         selection = questdlg(['Found version ' installedVersion ' of the Switched Mode Power Supply Toolbox installed. ' ...
        'The latest release from github is ' latestVersion '. Would you like to install the latest version?'], ...
        "Confirm Install");
        if strcmp(selection, 'Yes')
            matlab.addons.toolbox.uninstallToolbox(thisToolbox);
            installedToolbox = matlab.addons.toolbox.installToolbox(outfn);
        end
    end
else
    disp('Toolbox is up to date')
end

%% update component libraries
AURAdb(1).updateLibraries

