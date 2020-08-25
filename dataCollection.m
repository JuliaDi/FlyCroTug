% Created by Julia Di
% juliadi@stanford.edu
% modified August 18, 2020
%%
% Settings for  Arduino serial port object
% serialportlist("available")'
clear arduinoObj;
port = 'COM3';
board = 'Due';
baudRate = 115200;

% Create serial object for Arduino
arduinoObj = serial(port, 'BaudRate', baudRate, 'Terminator', 'LF');
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
% pixels per thou / 0.0254 thou per mm
dpm = 5.25 / 0.0254;

% Create figure for displacement vs force data
Smax = 45; % Total number of spines on cassette

% spine model
% stiffness is N / mm from data
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
plot(spines, sigmas,'or');
title('Spines vs. Sample Variance')
xlabel('Number of Spines'); ylabel('Stiffness Variance');
%stiff_interp = interp1(spines,stiffness,spinerange);
%plot(stiffness,spines,'o',stiff_interp, spinerange,':.');

%% plotting
Dmax = 0.5; % Total displacement for data collection (mm)
Fmax = 35; % max force for data collection (N)
x = linspace(0,Dmax, 200);
% live figure
figure,
grid on,
xlabel ('Displacement (mm)'), ylabel('Force Data'),
axis([0 Dmax 0 Fmax]),
hold on
plot(x, m_45*x, x, m_35*x, x, m_25*x, x, m_15*x),
% Read and plot the data from Arduino
Ts = 0.02; % Sampling time (s)
i = 0;
force = 0;
mouseCoord = [0, 0];
displacement = 0;
xdisplacement = 0;
sampled_stiffness = 0;
t = 0;
tic % Start timer
Tmax = 100; % Total time for data collection (s)


while toc <= Tmax
    i = i + 1; % because i starts at 0, and we will use i for indexing
    % Read arduino buffer data
    force(i) = abs(str2double(fscanf(arduinoObj)));
    %force(i) = fscanf(arduinoObj, '%e', 7);
    %force(i) = 1;
    
    % get current mouse position
    ans , mouse = get(0, 'PointerLocation'); % mouse is x,y coordinates
    mouseCoord(i+1,1) = mouse(1); mouseCoord(i+1,2) = mouse(2);
    % calculating displacement in the x direction
    xdisplacement(i) = abs(mouseCoord(i+1,1) - mouseCoord(2,1))./dpm;
    fprintf('x: '); disp(xdisplacement(i));
    fprintf('force: '); disp(force(i));
    %fprintf('y: '); disp(mouse(2));
    fprintf('-----\n');
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
        %plot([xdisplacement(i-1) xdisplacement(i)],[force(i-1) force(i)], '.');
        %drawnow;
        
        % calculate probability
        
%         % find the slope from samples
%         p = polyfit(xdisplacement(1:i),force(1:i),1); % get least squares line coefficients
%         sampled_stiffness(i) = p(1); % the sampled slope
%         k = dsearchn(stiff_interp', sampled_stiffness(i)); %  returns the indices of the closest points in stiff_interp to sampled stiffness in Euclidean distance.
%         pred_num_spines(i) = spinerange(k);
%         % add in confidence interval
%         fprintf('Predicted Number of Spines: ');
%         disp(pred_num_spines(i));
%         % fprintf('Confidence Interval: ');
%         % disp();
    end
end

%%
% close serial port when done
fclose(arduinoObj);
delete(arduinoObj);
clear arduinoObj;

%% Empty data for processing
concatenated_35 = [];
concatenated_35pixels = [];
%% Processing
force = force'; % make column vector
xdisplacement = xdisplacement'; % make column vector
k = 1;
pixel_trial = mouseCoord(2:end,1);
trial = horzcat(xdisplacement, force);
concatenated_35 = [concatenated_35 trial];
concatenated_35pixels = [concatenated_35pixels pixel_trial];

%%
%trial35 = force;
%pixel_trial35 = mouseCoord(2:end,1);
% pixel_trial35 = [pixel_trial35(1); pixel_trial35];
for i = 1:446
    disp_trial35(i) = abs(pixel_trial35(i+1) - pixel_trial35(1))./dpm;
end

disp_trial35 = disp_trial35';

%%
%force = force'; % make column vector
%trial35_2 = force;
%pixel_trial35_2 = mouseCoord(2:end,1);
%pixel_trial35_2 = [pixel_trial35_2(1); pixel_trial35_2];
%for i = 1:548
%    disp_trial35_2(i) = abs(pixel_trial35_2(i+1) - pixel_trial35_2(1))./dpm;
%end

disp_trial35_2 = disp_trial35_2';

%%
force = force'; % make column vector
xdisplacement = xdisplacement'; % make column vector
disp_trial35_3 = xdisplacement;
trial35_3 = force;
pixel_trial35_3 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial35_4 = xdisplacement;
trial35_4 = force;
pixel_trial35_4 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial25_1 = xdisplacement;
trial25_1 = force;
pixel_trial25_1 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial25_2 = xdisplacement;
trial25_2 = force;
pixel_trial25_2 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial25_3 = xdisplacement;
trial25_3 = force;
pixel_trial25_3 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial25_4 = xdisplacement;
trial25_4 = force;
pixel_trial25_4 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial30_1 = xdisplacement;
trial30_1 = force;
pixel_trial30_1 = mouseCoord(2:end,1);

%% displacement seemed off for this one
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial30_2 = xdisplacement;
trial30_2 = force;
pixel_trial30_2 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial30_3 = xdisplacement;
trial30_3 = force;
pixel_trial30_3 = mouseCoord(2:end,1);

%% force reading seemed low?
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial45_1 = xdisplacement;
trial45_1 = force;
pixel_trial45_1 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial45_2 = xdisplacement;
trial45_2 = force;
pixel_trial45_2 = mouseCoord(2:end,1);

%%
force = force';
xdisplacement = xdisplacement'; % make column vector
disp_trial45_3 = xdisplacement;
trial45_3 = force;
pixel_trial45_3 = mouseCoord(2:end,1);

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
%% 45 spines raw data
figure()
hold on
plot(disp_trial45_1(155:217), trial45_1(155:217), 'b');
b1 = [ones(217-155+1,1), disp_trial45_1(155:217)]\trial45_1(155:217); % linear regression 1
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b1, '--b');

plot(disp_trial45_2(116:195), trial45_2(116:195), 'r');
b2 = [ones(195-116+1,1), disp_trial45_2(116:195)]\trial45_2(116:195); % linear regression 2
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)'] * b2, '--r');

plot(disp_trial45_3(75:158), trial45_3(75:158), 'm');
b3 = [ones(158-75+1,1), disp_trial45_3(75:158)]\trial45_3(75:158); % linear regression 3
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)'] * b3, '--m');
plot(linspace(0,0.2,100)', m_45*linspace(0,0.2,100)', 'LineWidth',2);
legend('1','1 linear fit', '2','2 lin fit','3','3 lin fit', 'model slope','Location','southeast');
title('45 Spines raw data');
xlabel('Displacement (mm)'); ylabel('Force (N)');
hold off

%% 35 spines raw data
% plot force vs displacement for n = 35 spines
% force: concatenated_35(:,2), trial35, trial35_2, trial35_3, trial35_4
% displacement: concatenated_35(:,1), disp_trial35_3, disp_trial35_4

figure()
hold on
plot(concatenated_35(246:end,1), concatenated_35(246:end,2), 'r');
b1 = [ones(129,1), concatenated_35(246:end,1)]\concatenated_35(246:end,2); % linear regression 1
plot(linspace(0,0.5,100)', [ones(100,1), linspace(0,0.5,100)']*b1, '--r');

plot(disp_trial35(263:403), trial35(263:403)),'m';
b2 = [ones(403-263+1,1), disp_trial35(263:403)]\trial35(263:403);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b2, '--m');

plot(disp_trial35_2(200:472), trial35_2(200:472), 'b');
b3 = [ones(273,1), disp_trial35_2(200:472)]\trial35_2(200:472);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b3, '--b');

plot(disp_trial35_3(172:352), trial35_3(172:352),'k');
b4 = [ones(352-172+1,1), disp_trial35_3(172:352)]\trial35_3(172:352);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b4, '--k');

plot(disp_trial35_4(106:244), trial35_4(106:244), 'g');
b5 = [ones(244-106+1,1), disp_trial35_4(106:244)]\trial35_4(106:244);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b5, '--g');

plot(x, m_35*x, 'LineWidth',2);
legend('1','1 lin fit', '2', '2 lin fit', '3', '3 lin fit', '4', '4 lin fit', '5', '5 lin fit','model slope','Location','southeast');
title('35 Spines raw data');
xlabel('Displacement (mm)'); ylabel('Force (N)');
hold off
%% 30 spines raw data
figure()
hold on
plot(disp_trial30_1(292:403), trial30_1(292:403), 'r');
b1 = [ones(403-292+1,1), disp_trial30_1(292:403)]\trial30_1(292:403);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b1, '--r');

plot(disp_trial30_2(97:213), trial30_2(97:213),'m');
b2 = [ones(213-97+1, 1), disp_trial30_2(97:213)]\trial30_2(97:213);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b2, '--m');

plot(disp_trial30_3(122:188), trial30_3(122:188),'b');
b3 = [ones(188-122+1,1), disp_trial30_3(122:188)]\trial30_3(122:188);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b3, '--b');

plot(linspace(0,0.2,100)', stiff_interp(30+1)*linspace(0,0.2,100)', 'LineWidth', 2);
legend('1','1 lin fit', '2', '2 lin fit','3', '3 lin fit','interpolated model slope','Location','southeast');
title('30 Spines raw data');
xlabel('Displacement (mm)'); ylabel('Force (N)');
hold off

%% 25 spines raw data
figure()
hold on
plot(disp_trial25_1(641:777), trial25_1(641:777),'r');
b1 = [ones(777-641+1,1), disp_trial25_1(641:777)]\trial25_1(641:777);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b1, '--r');

plot(disp_trial25_2(138:297), trial25_2(138:297),'m');
b2 = [ones(297-138+1,1), disp_trial25_2(138:297)]\trial25_2(138:297);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b2, '--m');

plot(disp_trial25_3(119:224), trial25_3(119:224),'b');
b3 = [ones(224-119+1,1), disp_trial25_3(119:224)]\trial25_3(119:224);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b3, '--b');

plot(disp_trial25_4(106:263), trial25_4(106:263),'g');
b4 = [ones(263-106+1,1), disp_trial25_4(106:263)]\trial25_4(106:263);
plot(linspace(0,0.2,100)', [ones(100,1), linspace(0,0.2,100)']*b4, '--g');

plot(linspace(0,0.3,100)', m_25*linspace(0,0.3,100)', 'LineWidth', 2);
legend('1','1 lin fit', '2', '2 lin fit','3', '3 lin fit','4','4 lin fit','model slope','Location','southeast');
title('25 Spines raw data');
xlabel('Displacement (mm)'); ylabel('Force (N)');
hold off


%%
        % find the slope from samples
        p = polyfit(xdisplacement(1:end),force(1:end),1); % get least squares line coefficients
        sampled_stiffness = p(1); % the sampled slope
        k = dsearchn(stiff_interp', sampled_stiffness); %  returns the indices of the closest points in stiff_interp to sampled stiffness in Euclidean distance.
        pred_num_spines = spinerange(k);
        % add in confidence interval
        fprintf('Predicted Number of Spines: ');
        disp(pred_num_spines);
        % fprintf('Confidence Interval: ');
        % disp();

writematrix(M,'M.xls')