fMode = 'triggeredspec';

switch fMode
    case {'compute','triggeredspec'}
        display('hi')
end
switch fMode
    case 'triggeredspec'
        display('world')
end