function [km] = K(theta, k)
    U = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    Un = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    km=Un*k*U;
end
