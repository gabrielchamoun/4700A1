clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')

eCount = 1000;      % Total number of electrons
ePlotted = 10;      % Number of Electrons Plotted
dt = 10e-15;        % Time step 10fs -> ( -(Width/100) / vT )
tStop = 100 * dt;	% Stop Time

kB = 1.38066e-23;   % J/K 
m0 = 9.11e-31;
mn = 0.26*m0;

Width = 200e-9;
Height = 100e-9;

% Thermal Velocity
Temp = 300;     % K
vT = sqrt((2*kB*Temp)/mn);

% scattering
% probability of scattering
toggleS = 1;    % Toggle Scattering OFF(0) ON(1)
tmin = 0.2e-12;
pScatter = 1 - exp(-dt/tmin);

mfp = vT * tmin;
fprintf("tmin = %d m\n", tmin);
fprintf("Mean Free Path = %d m\n", mfp);


% bottlenecking
% Addition of box colliders to the sim
% Providing opposing corners of desired boxes
box = [
    0.8, 0.8, 1.2, 1.2;
    0, 0.4, 0.4, 0;
    ];

% Initializing all electrons with a position and velocity
eObj = struct('x', 0, 'y', 0, 'vx', 0, 'vy', 0, 'vm', 0);  % Electron Object Classification
eCol = rand(eCount,3);   % Random colours for plotting
for i = 1 : eCount
    eObj(i).x = rand()*Width;
    eObj(i).y = rand()*Height;
    eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vm = sqrt(eObj(i).vx^2 + eObj(i).vy^2);
end

t = 0;          % Init time
counter = 1;    % Init Counter
while t < tStop

    for i = 1 : eCount
        
        % updating position
        eObj(i).x(counter+1) = eObj(i).x(counter) + eObj(i).vx * dt;
        eObj(i).y(counter+1) = eObj(i).y(counter) + eObj(i).vy * dt;
        
        % scattering effect - Randomize direction/magnitude of velocity
        % probability of scattering based on p.
        if pScatter > rand() && toggleS	% 'if true'
            eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
            eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
        end
        % magnitude of velocity
        eObj(i).vm = sqrt(eObj(i).vx.^2 + eObj(i).vy.^2);
        
        % plotting only the first 10 electrons
        if (i <= 10)
            subplot(2,1,1)
            p = plot( [eObj(i).x(counter), eObj(i).x(counter+1)], ...
                [eObj(i).y(counter), eObj(i).y(counter+1)] );
            p.Color = eCol(i,:);
            hold on
        end
        
        % boundary conditions
        % y = 200nm boundary
        if eObj(i).x(counter+1) > Width
            eObj(i).x(counter+1) = eObj(i).x(counter+1) - Width;
        end
        
        % y = 0nm boundary
        if eObj(i).x(counter+1) < 0
            eObj(i).x(counter+1) = eObj(i).x(counter+1) + Width;
        end
        
        % x = 100nm boundary
        if eObj(i).y(counter+1) > Height
            diff = eObj(i).y(counter+1) - Height;
            eObj(i).y(counter+1) = Height - diff;
            eObj(i).vy = -eObj(i).vy;
        end
        
        % x = 0nm boundary
        if eObj(i).y(counter+1) < 0
            diff = -eObj(i).y(counter+1);
            eObj(i).y(counter+1) = diff;
            eObj(i).vy = -eObj(i).vy;
        end        
    end
    pause(0.05);       % Delay for animation
    axis([0,Width,0,Height]);  % Plot Axis' set
    t = t + dt;         % Incrementing Time
    
    Time(:,counter) = t;
    averageTemperature = ( ([eObj(:).vm].^2) .* mn ) ./ (kB*2);
    Temp(:,counter) = mean(averageTemperature); % Average Temperature
    subplot(2,2,3), plot(Time, Temp);
    subplot(2,2,4), histogram([eObj(:).vm],50)
    counter = counter + 1;      % Incrementing Sim Counter
end
hold off
