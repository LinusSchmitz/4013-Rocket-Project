function [Ab,Vb,Vc] = burn_geometry(r,h,rb, rout)

%  This burnback model only assumes burning with the cylindrical 
%  perforation of the grain.  Modify as needed for your rocket motor. 
 %outer burn radius

    if rb+r >= rout % motor is burnt out
        Ab = 0; % [m^2] 
    else % there is grain remaining
        %% BURN AREA
        %Ab = 2*pi*(r + rb)*h; % [m^2] total burn area for cylindrically perforated grain
        Ab_x = 2*(pi*rout^2 - pi*(r+rb)^2);
        Ab_y = 2*pi*(r+rb)*(h-2*rb);
        Ab = Ab_x + Ab_y;

        %% VOLUME of PROPELLANT CONSUMED
        Vb = pi*((r + rb)^2 - r^2)*(h-rb); % [m^3]

        %% CHAMBER VOLUME
        Vc = pi*(r + rb)^2*(h-rb) + 2*pi*rout^2*rb; % [m^3]  ignores the volume between the end of 
        % grain and nozzle entrance 
        
        end