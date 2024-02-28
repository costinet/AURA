gitRelease = webread('https://api.github.com/repos/costinet/aura/releases/latest');
latestVersion = gitRelease.tag_name;
installedVersion = matlab.addons.toolbox.toolboxVersion('switched')
if 
releaseURL = gitRelease.assets.browser_download_url;
outfn = websave(gitRelease.assets.name, releaseURL);
installedToolbox = matlab.addons.toolbox.installToolbox(outfn);