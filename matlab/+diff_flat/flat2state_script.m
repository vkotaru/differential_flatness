

flats.x = zeros(3,1);
flats.dx = zeros(3,1);
flats.d2x = zeros(3,1);
flats.d3x = zeros(3,1);
flats.d4x = zeros(3,1);

obj.m = 1;
obj.J = eye(3)*1e-2;
obj.e1 = [1; 0; 0];
obj.g = 9.81;
obj.e3 = [0; 0; 1];
obj.e2 = [0; 1; 0];

ref = flat2stateQuad(obj, flats);


