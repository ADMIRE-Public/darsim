function start_diary(folder)
diary off
diary_folder = fullfile(folder,'Output');
%% Set up Diary
if ~exist(diary_folder, 'dir')
    mkdir(diary_folder);
else
    rmdir(diary_folder, 's');
    mkdir(diary_folder);
end
diary_path = fullfile(diary_folder,"RunDiary.txt");
diary(diary_path);
end