clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')

eCount = 1000;      % Total number of electrons
ePlotted = 10;      % Number of Electrons Plotted
dt = 10e-15;        % Time step 10fs -> ( -(Width/100) / vT )
tStop = 200 * dt;	% Stop Time

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
fprintf("Nominal Mean Free Path = %d m\n", mfp);


% bottlenecking
% Addition of box colliders to the sim
% Providing opposing corners of desired boxes
box = [
    0.8, 0, 1.2, 0.4;
    0.8, 0.6, 1.2, 1;
    ].*1e-7;


% Initializing all electrons with a position and velocity
eObj = struct('x', 0, 'y', 0, 'vx', 0, 'vy', 0, 'vm', 0);  % Electron Object Classification
eCol = rand(eCount,3);   % Random colours for plotting
for i = 1 : eCount
    eObj(i).x = rand()*Width;
    eObj(i).y = rand()*Height;
    for b = 1:size(box,1)
        while eObj(i).x < Width && eObj(i).x > 0 && eObj(i).y < Height && eObj(i).y > 0 && ...
        eObj(i).x <= max(box(b,1),box(b,3)) && eObj(i).x >= min(box(b,1),box(b,3)) && ...
        eObj(i).y <= max(box(b,2),box(b,4)) && eObj(i).y >= min(box(b,2),box(b,4))
            eObj(i).x = rand()*Width;
            eObj(i).y = rand()*Height;
        end
    end
    eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vm = sqrt(eObj(i).vx^2 + eObj(i).vy^2);
end

t = 0;          % Init time
counter = 2;    % Init Counter
while t < tStop

    for i = 1 : eCount
        
        % updating position
        eObj(i).x(counter) = eObj(i).x(counter-1) + eObj(i).vx * dt;
        eObj(i).y(counter) = eObj(i).y(counter-1) + eObj(i).vy * dt;
        
        % scattering effect - Randomize direction/magnitude of velocity
        % probability of scattering based on p.
        if pScatter > rand() && toggleS	% 'if true'
            eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
            eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
        end
        
        
        % magnitude of velocity
        eObj(i).vm = sqrt(eObj(i).vx.^2 + eObj(i).vy.^2);
        
        % plotting only the first 10 electrons
        if (i <= ePlotted)
            subplot(3,1,1)
            p = plot( [eObj(i).x(counter-1), eObj(i).x(counter)], ...
                [eObj(i).y(counter-1), eObj(i).y(counter)] );
            p.Color = eCol(i,:);
            hold on
        end
        
        % boundary conditions
        if eObj(i).x(counter) > Width   % y = 200nm boundary
            eObj(i).x(counter) = eObj(i).x(counter) - Width;
        end
        
        if eObj(i).x(counter) < 0   % y = 0nm boundary
            eObj(i).x(counter) = eObj(i).x(counter) + Width;
        end
        
        if eObj(i).y(counter) > Height  % x = 100nm boundary
            diff = eObj(i).y(counter) - Height;
            eObj(i).y(counter) = Height - diff;
            eObj(i).vy = -eObj(i).vy;
        end
        
        if eObj(i).y(counter) < 0   % x = 0nm boundary
            diff = -eObj(i).y(counter);
            eObj(i).y(counter) = diff;
            eObj(i).vy = -eObj(i).vy;
        end
        
        for b = 1:size(box,1)
            if counter == 2 % Only plot once per simulation
                boxX = [box(b,1), box(b,1), box(b,3), box(b,3), box(b,1)];
                boxY = [box(b,2), box(b,4), box(b,4), box(b,2), box(b,2)];
                fill(boxX, boxY, 'k');
            end
            
            % Checking for collision with boundary
            if eObj(i).y(counter) <= max(box(b,2),box(b,4)) && ... 
            eObj(i).y(counter) >= min(box(b,2),box(b,4))
                if eObj(i).x(counter) <= max(box(b,1),box(b,3)) && ...
                eObj(i).x(counter) > max(box(b,1),box(b,3)) - (Width/10) || ...
                eObj(i).x(counter) >= min(box(b,1),box(b,3)) && ...
                eObj(i).x(counter) < min(box(b,1),box(b,3)) + (Width/10)
                    eObj(i).vx = -eObj(i).vx;
                end
            end
            if eObj(i).x(counter) <= max(box(b,1),box(b,3)) && ...
            eObj(i).x(counter) >= min(box(b,1),box(b,3))
                if eObj(i).y(counter) <= max(box(b,2),box(b,4)) && ...
                eObj(i).y(counter) > max(box(b,2),box(b,4)) - (Height/10) || ...
                eObj(i).y(counter) >= min(box(b,2),box(b,4)) && ...
                eObj(i).y(counter) < min(box(b,2),box(b,4)) + (Height/10)
                    eObj(i).vy = -eObj(i).vy;
                end
            end
        end
        
    end
    pause(0.05);       % Delay for animation
    axis([0,Width,0,Height]);  % Plot Axis' set
    t = t + dt;         % Incrementing Time
    
    % Plotting the Average Temperature over time
    Time(:,counter) = t;
    averageTemperature = ( ([eObj(:).vm].^2) .* mn ) ./ (kB*2);
    Temp(:,counter) = mean(averageTemperature); % Average Temperature
    subplot(3,1,2), plot(Time, Temp, 'k', 'LineWidth',1.75);
    
    % Plotting the distribution of velocities
    subplot(3,2,5), histogram([eObj(:).vm],50)
    
    % Plotting the MEAN FREE PATH over time
    mfp(:,counter) = mean([eObj(:).vm]) * tmin;
    subplot(3,2,6), plot(Time, mfp, 'k', 'LineWidth',1.75);
    counter = counter + 1;      % Incrementing Sim Counter
end
hold off
