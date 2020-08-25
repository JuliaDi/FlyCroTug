% Created by Julia Di
% juliadi@stanford.edu
% August 09, 2020
%%
% Settings for  Arduino serial port object
% serialportlist("available")'
clear arduinoObj;
port = 'COM4';
baudRate = 115200;

% Create serial object for Arduino
arduinoObj = serial(port, 'BaudRate', baudRate, 'Terminator', 'CR/LF');
try
    fopen(arduinoObj);
catch err
    fclose(instrfind);
    error('Check port configurations. Remember to clear when done.');
end
flushinput(arduinoObj);
%% Spine model

% calibration variables
% calibrate force on Arduino

% calibrate displacement
% 6.4 pixels per thou = 6.4 / 0.0254 pixels per mm
dpm = 6.4 / 0.0254;

% Create figure for displacement vs force data
Smax = 45; % Total number of spines on cassette

% spine model
% stiffness is N per mm from data
m_45 = 85.82;
m_35 = 82.46;
m_25 = 76.64;
m_15 = 68.54;

spines = [15, 25, 35, 45];
spinerange = 0 : Smax;
stiffness = [m_15, m_25, m_35, m_45];

% plot stiffness on x axis and number of spines on y axis
% then run linear interpolation to get intermediary points
figure()
stiff_interp = interp1(spines,stiffness,spinerange);
plot(stiffness,spines,'o',stiff_interp, spinerange,':.');
title('Linear Interpolation of Stiffness vs Spines');
xlabel('Stiffness (N/mm)'); ylabel('Number of Spines');

%stiff_interp gives the slope values for the intermediary # of spines
% Need to do linear interpolation to get the std dev at
% intermediary points
sigma_45 = 3.98;
sigma_35 = 3.75;
sigma_25 = 3.81;
sigma_15 = 4.96;
sigmas = [sigma_15, sigma_25, sigma_35, sigma_45];
figure()
plot(spines, sigmas,'or', 'Linewidth', 2);

sigma_interp = interp1(spines,sigmas,spinerange);
plot(spines,sigmas,'o',spinerange, sigma_interp,':.');
title('Training Spines vs. Sample Standard Deviation')
xlabel('Number of Spines'); ylabel('Stiffness Standard Deviation');


%% plotting
Dmax = 40; % Total displacement for data collection (mm)
x = 0:Dmax;
% live figure
figure,
grid on,
xlabel ('Displacement (mm)'), ylabel('Force Data'),
%xlabel ('Stiffness (N/mm)'), ylabel('Number of Spines'),
%axis([0 Dmax 0 17]),
hold on
% plot(x, y_45, x, y_35, x, y_25, x, y_15),
% Read and plot the data from Arduino
Ts = 0.02; % Sampling time (s)
i = 0;
force = 0;
mouseCoord = [0, 0];
displacement = 0;
sampled_stiffness = 0;
t = 0;
tic % Start timer
Tmax = 10; % Total time for data collection (s)

% for testing only
% 45 spines
%displacement = [0,0.0254,0.0508,0.0762,0.1016,0.127,0.1524,0.1778,0.2032,0.2286,0.254,0.2794,0.3048,0.3302,0.3556];
%force = [9.1,10.8,12.3,14.4,16.8,18.9,21.1,23.4,25.8,28,30.5,33,35.3,37.4,39.7];

% 35 spines
%displacement = [0,0.0254,0.0508,0.0762,0.1016,0.127,0.1524,0.1778,0.2032,0.2286,0.254,0.2794,0.3048,0.3302,0.3556,0.381,0.4064];
%force = [9.3,11.1,12.7,14.1,15.9,17.6,19.4,21.2,23.3,25.4,27.3,29.5,31.5,33.5,35.7,37.4,39.6];

% 15 spines
%displacement = [0,0.0254,0.0508,0.0762,0.1016,0.127,0.1524,0.1778,0.2032,0.2286,0.254,0.2794,0.3048,0.3302,0.3556,0.381];
%force = [8,9.7,11.2,13,14.7,16.2,18.2,20.1,22.3,24.3,26.3,28.3,30.4,32.1,34.1,36.1];


while toc <= Tmax
    i = i + 1; % because i starts at 0, and we will use i for indexing
    % Read arduino buffer data
    force(i) = fscanf(arduinoObj, '%f');
    
    % get current mouse position
    ans , mouse = get(0, 'PointerLocation'); % mouse is x,y coordinates
    mouseCoord(i+1,1) = mouse(1); mouseCoord(i+1,2) = mouse(2);
    % calculating total displacement
    displacement(i) = norm(mouseCoord(i+1,:)./dpm - mouseCoord(2,:)./dpm);
    %displacement(i) = sqrt(sum((mouseCoord(i+1,:)./dpm - mouseCoord(2,:)./dpm).^2)); %Euclidean distance from starting point
    
    % Read time stamp
    % If reading faster than sampling rate, force sampling time.
    % If reading slower than sampling rate, nothing can be done. Consider
    % decreasing the set sampling time Ts
    t(i) = toc;
    if i > 1
        T = toc - t(i-1);
        while T < Ts
            T = toc - t(i-1);
        end
    end
    
    t(i) = toc;
    
    % Plot live data
    if i > 1
        %line([t(i-1) t(i)],[data(i-1) data(i)]);
        plot([displacement(i-1) displacement(i)],[force(i-1) force(i)], '.');
        drawnow;
        
        % calculate probability
        
        % find the slope from samples
        p = polyfit(displacement(1:i),force(1:i),1); % get least squares line coefficients
        sampled_stiffness(i) = p(1); % the sampled slope
        
        % or find slope between datapoints?
        % sampled_stiffness(i) = (force(i) - force(i-1)) / (displacement(i) - displacement(i-1));
        
        
        k = dsearchn(stiff_interp', sampled_stiffness(i)); %  returns the indices of the closest points in stiff_interp to sampled stiffness in Euclidean distance.
        pred_num_spines(i) = spinerange(k);

        % add in confidence interval
        fprintf('Predicted Number of Spines: ');
        disp(pred_num_spines(i));
        fprintf('Confidence Interval: ');
        disp('90%');
        disp('------');
    end
end


% close serial port when done
fclose(arduinoObj);
delete(arduinoObj);
clear arduinoObj;




%% For MatLab 2019b
% %configureTerminator(arduinoObj, "CR/LF");
% flushinput(arduinoObj); % remove any old data
% fopen(arduinoObj)l
% fprintf(arduino, '%s', char(doi));
% fclose(arduino);
% arduinoObj.UserData = struct("Data",[],"Count", 1);
% 
% % window for detecting mouse coordinates
% figure(1)
% set (gcf, 'WindowButtonMotionFcn', @mouseMove);
% 
% % read force data
% configureCallback(arduinoObj, "terminator", @readForceData);
% clear;
% 
% % read mouse movements 
% function mouseMove (object, eventdata)
%     
%     C = get (gca, 'CurrentPoint');
%     x = C(1,1); y = C(1,2);
%     title(gca, ['(X,Y) = (', num2str(x), ', ',num2str(y), ')']);
%     
% end
% 
% function readForceData(src, ~)
% 
% data = readline(src);
% 
% % convert string data to double
% src.UserData.Data(end+1) = str2double(data); 
% 
% src.UserData.Count = src.UserData.Count + 1;
% 
% %figure(1)
% %plot(src.UserData.Data(2:end));
% %drawnow;
% 
% if src.UserData.Count > 500
%     configureCallback(src, "off");
% end
% end
% 
% 
% 
