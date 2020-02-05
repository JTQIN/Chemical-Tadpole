function clocktime = clocktime2text
%
% KLF: Turn the clocktime into text.  Useful for making unique filenames
% when saving files.

temptime = clock;    % Vector [year month day hour minute second]
tempname = [num2str(temptime(1)) num2str(temptime(2),'%02d') ...
    num2str(temptime(3),'%02d'), num2str(temptime(4),'%02d') ...
    num2str(temptime(5),'%02d'), num2str(round(temptime(6)),'%02d')];
clocktime = num2str(tempname);