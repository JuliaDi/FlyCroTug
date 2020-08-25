clear a;
try
    a = arduino('COM4', 'Mega2560', 'Libraries', 'Servo');
catch err
    clear a;
    error('Check port configurations. Remember to clear when done.');
end

%%

% from datasheet:
% min pulse: 500 microseconds
% max pulse: 2500 microseconds
% https://github.com/microrobotics/DS3235-270/blob/master/DS3235-270_datasheet.pdf
clear s;
try 
    s = servo(a, 'D2', 'MinPulseDuration', 500*10^-6, 'MaxPulseDuration', 2500*10^-6);
catch err 
    clear s;
    error('Check servo configurations.');
end

%%

% read and write servo angle
% servo goes from 0 to 1 duty cycle -> 0 to 270 degrees
for angle = 0:0.2:1
    writePosition(s, angle);
    current_pos = readPosition(s);
    current_pos = current_pos*270;
    fprintf('Current motor position is %d degrees\n', current_pos);
    pause(2);
end

%%
clear s a