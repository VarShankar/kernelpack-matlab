function paths = download_rbc_capstone_data(source)
%DOWNLOAD_RBC_CAPSTONE_DATA Download and verify the published RBC dataset.

arguments
    source (1, 1) string = "https://github.com/VarShankar/kernelpack-matlab/releases/download/data-v1/kernelpack-matlab-rbc-data-v1.zip"
end

root = fileparts(mfilename('fullpath'));
dataRoot = fullfile(root, 'data', 'ibamr_rbc_3d');
trajectoryPath = fullfile(dataRoot, 'trajectory');
resultMatPath = fullfile(root, 'moving_surface_adr_tp_rbc_capstone.mat');
resultCsvPath = fullfile(root, 'moving_surface_adr_tp_rbc_capstone.csv');

paths = struct('trajectory', trajectoryPath, ...
    'resultMat', resultMatPath, 'resultCsv', resultCsvPath);
if isfolder(trajectoryPath) && isfile(resultMatPath) && isfile(resultCsvPath)
    return;
end

if isfolder(trajectoryPath) || isfile(resultMatPath) || isfile(resultCsvPath)
    error('kp:data:PartialRBCDataset', ...
        ['A partial RBC dataset already exists. Remove the incomplete files ' ...
         'before downloading the release asset again.']);
end

expectedHash = "1fe40d02155b2ab9cb041211cbf4929d40ce6b382f3bd56fcdcb575ab8466344";
archivePath = [tempname, '.zip'];
extractPath = tempname;
mkdir(extractPath);
cleanup = onCleanup(@() cleanupTemporaryFiles(archivePath, extractPath));

fprintf('Downloading RBC capstone data...\n');
if isfile(source)
    copyfile(source, archivePath);
else
    websave(archivePath, source);
end
actualHash = sha256File(archivePath);
if actualHash ~= expectedHash
    error('kp:data:RBCChecksumMismatch', ...
        'RBC data checksum mismatch: expected %s, received %s.', ...
        expectedHash, actualHash);
end

unzip(archivePath, extractPath);
if ~isfolder(dataRoot)
    mkdir(dataRoot);
end
movefile(fullfile(extractPath, 'trajectory'), trajectoryPath);
movefile(fullfile(extractPath, 'moving_surface_adr_tp_rbc_capstone.mat'), ...
    resultMatPath);
movefile(fullfile(extractPath, 'moving_surface_adr_tp_rbc_capstone.csv'), ...
    resultCsvPath);
fprintf('RBC capstone data installed in %s\n', dataRoot);
end

function hash = sha256File(path)
fileId = fopen(path, 'rb');
if fileId < 0
    error('kp:data:RBCArchiveReadFailure', 'Unable to read %s.', path);
end
cleanup = onCleanup(@() fclose(fileId));
digest = java.security.MessageDigest.getInstance('SHA-256');
while true
    block = fread(fileId, 1024 * 1024, '*uint8');
    if isempty(block)
        break;
    end
    digest.update(typecast(block, 'int8'));
end
bytes = typecast(digest.digest(), 'uint8');
hash = lower(string(reshape(dec2hex(bytes, 2).', 1, [])));
end

function cleanupTemporaryFiles(archivePath, extractPath)
if isfile(archivePath)
    delete(archivePath);
end
if isfolder(extractPath)
    rmdir(extractPath, 's');
end
end
