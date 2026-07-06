function startup_darsim(varargin)
addpath(genpath(pwd))
persistent flag_init
if ~isempty(varargin)
    flag_init = varargin{1};
end
if isempty(varargin)
    if isempty(flag_init)
        flag_init = 1;
        clc
        displayFile("welcome.md")
        displayFile("authors.md")
        displayFile("repository.md")
    elseif flag_init == 1
        clc
        displayFile("resume.md");
    else
        flag_init = 1;
    end
end
end