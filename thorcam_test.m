
str = 'C:\Program Files\Thorlabs\Scientific Imaging\Scientific Camera Support\Scientific Camera Interfaces\SDK\DotNet Toolkit\dlls\Managed_64_lib';


% Change directory to where stored the dlls
cd(str);
str=[pwd filesep 'Thorlabs.TSI.TLCamera.dll'];

% Check is the DLL's exist
if ~exist(str,'file')
    warning('Could not located the DLL to run the ThorCam! Exiting software.');
    return
end

% Load the NET assembly framework
try
    fprintf('Loading the NET assembly ... ');
    asmInfo=NET.addAssembly(str);
    asmInfo.Classes;
    disp('Dot NET assembly loaded.');
catch ME
    warning('unable to load NET assemy')
    return;
end

% Open the SDK
try
    fprintf('Opening the camera SDK...');
    tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;
    disp('loaded.');
catch ME
    warning('unable to load camera sdk')
end

try
    fprintf('looking for cameras');    
    serialNumbers = tlCameraSDK.DiscoverAvailableCameras;
    disp('cameras found');
catch ME
    disp(ME)
    warning('unable to find cameras')
    tlCameraSDK.Dispose;                    % Unload the SDK
    return;
end

tlCameraSDK.Dispose;                    % Unload the SDK

disp('done')
