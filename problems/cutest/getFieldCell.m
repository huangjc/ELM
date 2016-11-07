function out = getFieldCell( record, convert, varargin )
% convert: toggle convertion cell2mat if it is 1

if isempty( convert )
    convert = 1;
end

out = cell(length(varargin),1);
for i = 1:length(varargin)
    fieldname = strsplit( varargin{i}, '.' );

    [np,nm] = size( record );

    T = cell(np,nm);
    for kprob = 1:np
        for kmthd = 1:nm
            rcrd = record{kprob, kmthd};

            if ( length(fieldname) == 1 )
                T{kprob,kmthd} = rcrd.(fieldname{1});
            elseif ( length(fieldname) == 2 )
                T{kprob,kmthd} = rcrd.(fieldname{1}).(fieldname{2});
            end
        end
    end
    
    if convert == 1
        % try to tranform T into matrix if elements has the same data type
        T = cell2mat(T);
    end
        
    if length(varargin) == 1
        out = T;
    else
        out{i} = T;
    end
end

