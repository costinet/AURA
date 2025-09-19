function [ltspicePath, ltspiceEXE]  = findLTspiceInstallDirectory(obj)
% findLTspice Attempts to find the LTspice installation directory.
%   [ltspicePath, ltspiceEXE]  = findLTspiceInstallDirectory(obj) returns the path to the LTspice 
%   installation directory if found, otherwise returns an empty string.
 
ltspicePath = '';
ltspiceEXE = '';
 
if ispc
    % --- Windows ---
    % Common installation directories for different versions
    progFiles = getenv('ProgramFiles');
    localAppData = getenv('LocalAppData');
    
    possiblePaths = {
        fullfile(progFiles, 'ADI', 'LTspice'); ...
        fullfile(localAppData, 'Programs', 'ADI', 'LTspice'); ...
        fullfile(progFiles, 'LTC', 'LTspiceXVII'); ...
        fullfile(getenv('USERPROFILE'), 'Documents', 'LTspiceXVII'); ...
    };
    
    % Check for the existence of the main executable
    for i = 1:length(possiblePaths)
        if isfile(fullfile(possiblePaths{i}, 'XVIIx64.exe'))
            ltspicePath = possiblePaths{i};
            ltspiceEXE = fullfile(ltspicePath, 'XVIIx64.exe');
            return
        elseif isfile(fullfile(possiblePaths{i}, 'LTspice.exe'))
            ltspicePath = possiblePaths{i};
            ltspiceEXE = fullfile(ltspicePath, 'LTspice.exe');
            return;
        end
    end
    
    % Fallback: search the registry for uninstaller information
    % [~, result] = system('reg query "HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall" /s /f "LTspice"');
    % if ~isempty(result)
    %     matches = regexp(result, 'InstallLocation\s+REG_SZ\s+(.*)', 'tokens');
    %     if ~isempty(matches) && ~isempty(matches{1})
    %         ltspicePath = strtrim(matches{1}{1});
    %         if isfolder(ltspicePath)
    %              if isfile(fullfile(ltspicePath, 'XVIIx64.exe'))
    %                 ltspiceEXE = fullfile(ltspicePath, 'LTspice.exe');
    %                 return
    %             elseif isfile(fullfile(ltspicePath, 'LTspice.exe'))
    %                 ltspiceEXE = fullfile(ltspicePath, 'LTspice.exe');
    %                 return;
    %             end
    %             return;
    %         end
    %     end
    % end
 
elseif ismac
    % --- macOS ---
    % Default application location
    if isfolder('/Applications/LTspice.app')
        ltspicePath = '/Applications/LTspice.app';
        return;
    end
    
    % Check user's Applications folder
    userAppDir = fullfile(getenv('HOME'), 'Applications', 'LTspice.app');
    if isfolder(userAppDir)
        ltspicePath = userAppDir;
        return;
    end
    
    % % Check for user-specific library files which might indicate installation
    % userLibPath = fullfile(getenv('HOME'), 'Library', 'Application Support', 'LTspice');
    % if isfolder(userLibPath)
    %     % This is not the installation path, but confirms it's likely installed.
    %     % The primary executable is still the .app bundle. A more robust search
    %     % could be initiated from here if needed.
    %     [~, result] = system('mdfind "kMDItemFSName == ''LTspice.app''"');
    %     if ~isempty(result)
    %         ltspicePath = strtrim(result);
    %         return;
    %     end
    % end
 
elseif isunix
    % --- Linux (with Wine) ---
    winePath = fullfile(getenv('HOME'), '.wine', 'drive_c');
    
    if isfolder(winePath)
        % Common installation paths within Wine
        possibleWinePaths = {
            fullfile(winePath, 'Program Files', 'LTC', 'LTspiceXVII'); ...
            fullfile(winePath, 'Program Files', 'ADI', 'LTspice'); ...
        };
        
        for i = 1:length(possibleWinePaths)
            if isfile(fullfile(possibleWinePaths{i}, 'XVIIx64.exe'))
                ltspicePath = possibleWinePaths{i};
                ltspiceEXE = fullfile(ltspicePath, 'XVIIx64.exe');
                return
            elseif isfile(fullfile(possibleWinePaths{i}, 'LTspice.exe'))
                ltspicePath = possibleWinePaths{i};
                ltspiceEXE = fullfile(ltspicePath, 'LTspice.exe');
                return;
            end
        end
        
        % Fallback: search for the executable within the wine C: drive
        % [~, result] = system(['find "', winePath, '" -name "XVIIx64.exe" -o -name "LTspice.exe" 2>/dev/null | head -n 1']);
        % if ~isempty(result)
        %     ltspicePath = fileparts(strtrim(result));
        %     return;
        % end
    end
end
 
end
 