function out = biodex_load_file(filein)

if nargin < 1
    [filename, pathname] = uigetfile('*.dat', 'Data file...');
    filein = [pathname filename];
end

[header, data] = hdrload(filein);

field_names = {'MODE';'Away Limit (deg)';'Toward Limit (deg)';'Anatomic Ref (deg)';'Speed Away (deg/s)';...
    'Speed Toward (deg/s)';'Velocity Scaling (deg/s)';'Torque Scaling (Nm)';'Position Scaling';...
    'A_D Rate (Hz)'};

for i = 1:length(field_names)
    clear D
    fname = field_names{i};
    for j = 1:size(header,1)
        a = strfind(header(j,:),fname);
        if ~isempty(a)
            b = find(diff(isspace(header(j,:)))>0);
            D = header(j,length(fname)+2:b(end));
            % find any spaces
            e = findstr(fname, ' ');
            if ~isempty(e)
                fname(e) = '_';
            end
            % find any -
            e = findstr(fname, '-');
            if ~isempty(e)
                fname(e) = [];
            end
            %find any (
            e = findstr(fname, '(');
            if ~isempty(e)
                fname(e(1)-1:end) = [];
            end
            %find any /
            e = findstr(fname, '/');
            if ~isempty(e)
                fname(e) = [];
            end

            %find any ,
            e = findstr(fname, ',');
            if ~isempty(e)
                fname(e) = [];
            end

            %find any #
            e = findstr(fname, '#');
            if ~isempty(e)
                fname(e) = [];
            end

            %find any #
            e = findstr(fname, '%');
            if ~isempty(e)
                fname(e) = 'P';
            end
            if isnan(str2double(D))
                E = D;
            else E = str2double(D);
            end
            out.(fname) = E;
        end
    end
end

out.data = data;
    