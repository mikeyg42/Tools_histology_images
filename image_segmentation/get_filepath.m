function filepath = get_filepath(directory, basename)
    files = dir(fullfile(directory));
    filenames = {files.name};
    matches = regexp(filenames, ['^', basename, '\.\w+$'], 'match');
    match_idx = ~cellfun('isempty', matches);
    if any(match_idx)
        full_filename = matches{match_idx};
        filepath = fullfile(directory, full_filename);
       % [pathstr, filename, extension] = fileparts(full_filename);
    else
        error(['No file found matching basename ', basename]);
    end
end
