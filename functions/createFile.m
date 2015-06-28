function[] = createFile(f)
if (exist(f, 'file') == 2)
    delete(f); % delete file if it exists
end
diary(f); % create file to write command line output
end