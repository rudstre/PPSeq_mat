function bool = ishex(s)
    for i = 1:length(s)
        bool(i) = ~isempty(regexpi(s(i), '^#*[0-9A-Fa-f]+$', 'once'));
    end
end
