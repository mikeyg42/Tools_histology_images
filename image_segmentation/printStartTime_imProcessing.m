function printStartTime_imProcessing(activeFilename)
% prints in the command window the current date and time.

timeStr1 = datestr(now,'mm/dd');
disp(strcat('Today is : ',timeStr1));
timeStr2 = datestr(now,  'HH:MM:SS' );
disp(strcat('Processing for : ', activeFilename,' began at : ',timeStr2));

% adapted from the work by Yair of the Undocumented Matlab website
end
