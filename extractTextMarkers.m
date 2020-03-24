function stEntries = extractTextMarkers(sFilename)

    fid = fopen(sFilename);

    tmp = textscan(fid , '%s');
    stEntries = struct('start', [], 'end', [], 'source', []);

    for iEntry = 1:3:length(tmp{1})

        stEntries(end+1).start = str2double(tmp{1}{iEntry});
        stEntries(end).end = str2double(tmp{1}{iEntry+1});
        stEntries(end).source = tmp{1}{iEntry+2};

    end

    stEntries(1) = [];

    fclose(fid);

end

