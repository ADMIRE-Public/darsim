function darsim_menu()
menu_text = {'Is it your first time at Darsim? ';
    'Do you want to run a simulation? '};
menu_command = {['Learn with Tutorials'], ['learn_menu']; ...
    ['Find a case'], ['startup_darsim(2); DARSim(pwd);']; ...
    };
processed_command = build_menu_command(menu_command);
length_menu = length(menu_command);
display_bar
for index = 1:length_menu
    current_line = join([menu_text{index} processed_command{index}]);
    disp(current_line)
end
display_bar
end


function commands = build_menu_command(menu_command)
num_commands = size(menu_command,1);
commands = cell(size(menu_command,1),1);
for index = 1:num_commands
    commands{index} = create_link_menu(menu_command{index,1}, ...
                menu_command{index,2});
end
end


function string_link = create_link_menu(text, command, varargin)
if isempty(varargin)
    flag = false;
else
    flag = varargin{1};
end
if ~flag
    string_link = [create_matlablink(text,command)];
else
    string_link = [create_weblink(text,command)];
end
end