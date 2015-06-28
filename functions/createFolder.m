function[] = createFolder(outputFolderName)
if (exist(outputFolderName, 'dir'))
    rmdir(outputFolderName,'s');
end
    [status,~,~] = mkdir(outputFolderName);
    assert(status == 1);
end
