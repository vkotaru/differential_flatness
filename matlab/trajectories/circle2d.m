function [traj] = circle2d(t, varargin)

r = varargin{1}; % radius
w = varargin{2}; % angular velocity
c = varargin{3}; % center

traj.x = c + r*[sin(w*t); cos(w*t); 0];
traj.dx = r*[w*cos(w*t); -w*sin(w*t); 0];
traj.d2x = r*[-w^2*sin(w*t); -w^2*cos(w*t); 0];
traj.d3x = r*[-w^3*cos(w*t); w^3*sin(w*t); 0];
traj.d4x = r*[w^4*sin(w*t); w^4*cos(w*t); 0];

end