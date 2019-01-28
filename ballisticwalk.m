function ballisticwalk
% Perform forward simulation of ballistic walking  
  
% Instructions:
%
% 1. Use the Dynamics Workbench to generate equations of motion and
%    export them to an m-file "fballwalk" for state-derivative function.
%
% 2. Copy fballwalk into the present file, as a nested function. This
%    makes the parameter values (defined below) accessible to fballwalk.
%
% 3. Write an event function "eventballwalk" to detect the knee
%    reaching full extension, and use it to terminate the ode45 integration.
%
% 4. Copy the expressions for kinetic and potential energy into a 
%    separate function energyballwalk.
%
% 5. Run a forward simulation and verify energy conservation.


% some initial states, in order q1-q3, u1-u3   
x0 = [14; -14; -60; -50; 250; -150] * pi/180; 
% above we convert deg or deg/s into rad or rad/s
   
% Set the parameter values for simulation
l1 = 0.5; l2 = 0.5; lc1 = 0.433*l1; lc2 = 0.437*l2; % segment lengths, thigh and shank
m1 = 0.097; m2 = 0.06; % segment masses as fraction of body mass
BodyMass = 65; % convert to actual mass
M0 = (m1 + m2) * BodyMass; M1 = m1 * BodyMass; M2 = m2 * BodyMass;
eta1 = 0.54*l1; eta2 = 0.735*l2; % locations of center of mass
l0 = l1 + l2; lc0 = (m1*lc1 + m2*(l1+lc2)) / (m1 + m2); % stance leg length
I0 = M1*eta1^2 + (M2*eta2^2 + M2*l1^2) - (M1+M2)*lc0^2; % moments of inertia
I1 = M1 * eta1^2 - M1*lc1^2; I2 = M2 * eta2^2 - M2*lc2^2;
g = 9.81; % gravity  

% simulate ballistic walking forward for a short duration
options = odeset('events', @eventballwalk);
[ts, xs] = ode45(@fballwalk, [0 0.5], x0);

% Here is a plot of the segment angles, defined ccw from vertical
clf; subplot(131)
plot(ts, xs(:,1:3)); 
xlabel('time'); ylabel('segment angle');
legend('q1', 'q2', 'q3');

% Now check energy conservation, as a post-processing step
energies = zeros(length(ts),1);
kineticEnergies = zeros(length(ts),1);
gravPotentialEnergies = zeros(length(ts),1);

for i = 1:length(ts)
  [energies(i), kineticEnergies(i), gravPotentialEnergies(i)] = ...
    energyballwalk(xs(i,:));
end

% show all three energies
subplot(132)
plot(ts, energies, ts, kineticEnergies, ts, gravPotentialEnergies);
xlabel('time'); ylabel('energy');
legend('total energy', 'kinetic energy', 'potential energy');

% and zoom in on just total energy
subplot(133)
plot(ts, energies);
xlabel('time'); ylabel('energy');
legend('total energy');

% END OF MAIN FUNCTION (NESTED FILES BELOW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fballwalk goes below

% PUT CODE FOR FBALLWALK HERE

end % function

%% fballwalk should go above here, nested inside ballisticwalk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% energyballwalk computes kinetic and potential energy
function [energy, kineticEnergy, gravPotentialEnergy] = energyballwalk(x)
  
% COPY ENERGY EXPRESSIONS HERE FROM FBALLWALK
% BE SURE TO ALSO INCLUDE STATE ASSIGNMENTS AND TRIGONOMETRIC ASSIGNMENTS  

% State assignments

% Trig assignments (shortcuts for equations)

kineticEnergy = % CODE HERE

gravPotentialEnergy = % CODE HERE

energy = kineticEnergy + gravPotentialEnergy;

end % energyballwalk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, direction, isterminal] = eventballwalk(x)

% WRITE A FUNCTION TO DETECT KNEE FULL EXTENSION

% State assignments
q1 = x(1); q2 = x(2); q3 = x(3); 
u1 = x(4); u2 = x(5); u3 = x(6); 

% CODE HERE

end % eventballwalk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % ballisticwalk file
