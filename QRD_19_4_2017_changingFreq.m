%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Joey Hook - Special Project - 19/4/2017
%   Modelling the QRD
%   The Quadratic Residue Diffusor is famous!
%   This program makes a 2D version of one, then bounces sound
%   off and measures the reflection.
%   -Input signal is variable by user
%   -Effective diffusion bandwidth is variable by user
%       these parameters determine well width and depth
%   -Length of sequence and # of repetitions variable by user
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear, close
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input parameters:
%   f_lo = design frequency - lowest frequency efffectively diffused
%   f_hi = highest frequency effectively diffused
%   nW = number of wells in sequence
%   S  = number of times sequence is repeated
%   c  = wave speed (speed of sound in m/s)
%   L_ = length of each side (meters)
%   Fs = sampling rate
%   T  = total number of samples in simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f_lo = 400;
f_hi = 1000;
nW = 7;
S  = 6;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% flags
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   input: 'Sine', 'Nois', 'Hann', 'Gaus', 'Delt'
input = 'Sine'; 
% 0 = no picture, 1 = picture
animate = 1;

% run the program multiple times with increasing input frequency, 
% then plot all of them together
repeat = 9;
freqIN = zeros(1,repeat);
step = (f_hi-f_lo)/6;
for n=1:repeat
    freqIN(n) = f_lo+(n-2)*step;
end

Lx = 15;
Ly = Lx/2;

c  = 343;
SR = 15000;%f_hi*12*sqrt(2)
T  = 1000;

if SR>30000
    warning('SR is quite large, may run pretty slow')
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Derived parameters:
%   k      = seconds per timestep
%   lambda = courant # (unitless)
%   h      = width of gridspace (meters)
%   N_     = # of gridspaces in _ direction
%   w      = width of well (in gridspaces)
%   d      = length of divider (max depth of well in gridspaces)
%   wD     = total width of diffuser (in gridspaces)
%   l_sq   = courant # squared
%   Tf     = final time (in seconds)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
k = 1/SR;
lambda = 1/sqrt(2);
h = c*k/lambda;
Nx = floor(Lx/h);
Ny = floor(Ly/h);
w  = ceil(c/(2*f_hi)/h)
if w < 5
    warning('The fins are too large compared to well width.')
end
d  = ceil(c/(2*f_lo)/h);
wD = w*nW*S;
if wD > Nx
    error('The diffuser is bigger than the room.')
end
h = Lx/Nx % correction
lambda = c*k/h; % correction
l_sq = lambda^2;

% print diffuser dimensions in meters
width = w*h
depth = d*h
totalwidth = wD*h

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Build Diffuser :
% div  = position along x of divider (one fewer than # of wells)
% bot  = depth of wells
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
div = zeros(S*nW,1);
% equally spaced dividers
div(1:S*nW,1) = floor((0:S*nW-1)*w+0.5*Nx-0.5*wD);

bot = zeros(S*nW,1);
% quadratic residue sequence determines bottom of well
bot(1:S*nW) = mod(((1:S*nW)-1).^2,nW);

% make well depth same as largest number in sequence
mult = round(d/max(bot));
d = mult*max(bot);
bot = d-mult*bot;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initial conditions
%   u2 = 2 time steps ago
%   u1 = 1 time step ago
%   u  = present
%   P  = y position of input signal
%   Q  = x position of input signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% planes initialized at 3 time steps
u2 = zeros(Ny+1,Nx+1);
u1 = u2;
u  = u1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Boundry Conditions
% Matrix J: true=air, false=boundry
% Matrix K: Weighting based on J used in wave equation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
J = ones(Ny+1,Nx+1);
% outer walls
J(1,:) = 0; J(Ny+1,:) = 0; J(:,1) = 0; J(:,Nx+1) = 0;

% left and right most ends of diffuser
b1 = div(1)-round(w/2);
b2 = div(S*nW)+round(w/2);

% dividers
J(2:d+1,div(1:S*nW)) = 0;
% well bottoms
for n=2:S*nW
    J(bot(n)+1,div(n-1):div(n)) = 0;
end

% lid on diffusor
J(d+1,b1:b2) = 0;

% weighting matrix for wave equation
K = zeros(Ny+1,Nx+1);
K(2:Ny,2:Nx) = 2-l_sq*(J(3:Ny+1,2:Nx)+J(1:Ny-1,2:Nx)...
    +J(2:Ny,3:Nx+1)+J(2:Ny,1:Nx-1));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input/Output
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% position input signal in line with center of diffuser at opposite wall
P = round(Ny/2);
if mod(length(div),2)==1
    Q = div(ceil(length(div)/2));
else
    Q = div(ceil(length(div)/2))+ceil(w/2);
end

% output array radius
r = round(Ny-d);

% duration of input signal
f_Tf = (r-d)*h/c/3;
f_Nf = round(f_Tf/k);
time = 0:k:f_Tf;
% input signals f
switch (input)
    case {'Gaus'}
        bt = 0.3;
        span = 4;
        sps = 20;
        f = gaussdesign(bt,span,sps);
        repeat = 1;
    case {'Delt'}
        f = zeros(f_Nf,1);
        f(1)=1;
        repeat = 1;
    case {'Hann'}
        f(1:f_Nf) = 0.5*(1-cos(2*pi*(1:f_Nf)/f_Nf));
        repeat = 1;
    case {'Sine'}
        f = sin(2*pi*freqIN(1)*time);
        repeat = 9;
    otherwise
        f = 2*rand(1,f_Nf)-1;
        repeat = 1;
end
% fade out end of input
winlen = round(f_Tf/k/20);
halfHann(1:winlen) = -0.5*(1-cos(pi*(1:winlen)/winlen))+1;
f((length(f)+1-winlen):length(f)) = halfHann(1:winlen).*f((length(f)+1-winlen):length(f));

% input function factor
khsq = (k/h)^2;
maxkhsqft = max(khsq*f);

% samples for wave to travel across room and bounce off diffuser
bounce = (P)*h/c/k+length(f);

% initialize outputs
theta = 0:pi/36:pi;
numouts = length(theta);
outputsflat = zeros(T,numouts);
maxoutflat = zeros(numouts,repeat);
energyflat = zeros(numouts,repeat);
outputs = zeros(T,numouts);
maxout = zeros(numouts,repeat);
energy = zeros(numouts,repeat);
bigger=zeros(repeat,1);

% total duration
T = (2*r)*h/c/k+length(f);
Tf = T/SR;

% factor in wave eqn
fac = 1/(1+lambda);
for pass=1:repeat
    
% lid on diffusor
J(d+1,b1:b2) = 0;

% weighting matrix for wave equation
K = zeros(Ny+1,Nx+1);
K(2:Ny,2:Nx) = 2-l_sq*(J(3:Ny+1,2:Nx)+J(1:Ny-1,2:Nx)...
    +J(2:Ny,3:Nx+1)+J(2:Ny,1:Nx-1));

%update frequency
if input == 'Sine'
    f = sin(2*pi*freqIN(pass)*time);
    f((length(f)+1-winlen):length(f)) = halfHann(1:winlen).*f((length(f)+1-winlen):length(f));
end

%initialize grids
u2 = zeros(Ny+1,Nx+1);
u1 = u2;
u  = u1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% main loop - with flat plate
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tstart = tic;
for t=1:T
    
    % 2D wave equation
    u(2:Ny,2:Nx) = l_sq *( u1(3:Ny+1,2:Nx) + u1(1:Ny-1,2:Nx) ...
        + u1(2:Ny,3:Nx+1) + u1(2:Ny,1:Nx-1))...
        + u1(2:Ny,2:Nx).*(K(2:Ny,2:Nx))     ...
        - u2(2:Ny,2:Nx);
    
    % absorb wave at bottom wall either side of diffuser
    u(2,2:div(1)) = fac*(2*(1-l_sq)*u1(2,2:div(1)) ...
        + 2*l_sq*(u1(3,2:div(1)))...
        + lambda*u2(2,2:div(1)) ...
        - u2(2,2:div(1)));
    u(2,div(S*nW):Nx) = fac*(2*(1-l_sq)*u1(2,div(S*nW):Nx) ...
        + 2*l_sq*(u1(3,div(S*nW):Nx))...
        + lambda*u2(2,div(S*nW):Nx) ...
        - u2(2,div(S*nW):Nx));
    
    % absorb wave at top wall
    u(Ny,:) = fac*(2*(1-l_sq)*u1(Ny,:) + 2*l_sq*(u1(Ny-1,:))...
        + lambda*u2(Ny,:) - u2(Ny,:));
    
    % wait to turn on absorption on the side walls, and microphone array
    % until the input has bounced off the diffusor
    if t > bounce
        % absorb wave at left wall
        u(:,2) = fac*(2*(1-l_sq)*u1(:,2) + 2*l_sq*(u1(:,3))...
            + lambda*u2(:,2) - u2(:,2));
        
        % absorb wave at right wall
        u(:,Nx) = fac*(2*(1-l_sq)*u1(:,Nx) + 2*l_sq*(u1(:,Nx-1))...
            + lambda*u2(:,Nx) - u2(:,Nx));
           
        % output
        for thetaOUT=1:numouts
            outputsflat(t,thetaOUT) = u(abs(round(r*sin(theta(thetaOUT))+d)),round(r*cos(theta(thetaOUT))+Q));
            % here you can run signal in through your output locations if
            % you like.
            %u(abs(round(r*sin(theta(thetaOUT))+d)),round(r*cos(theta(thetaOUT))+Q))=khsq*(f(t)); 
        end
    end
    
    % zero value at boundries
    for n=b1:b2
        for m=2:d+1
            if J(m,n)==0
                u(m,n)=0;
            end
        end
    end
    
    % input function
    if t <= length(f)
        u(P,:) = khsq*f(t);
    end
    
    % advance previous two time steps
    u2 = u1;
    u1 = u;
    
    % draw
    if t==1 && pass==1 && animate==1
        scrsz = get(0,'ScreenSize');
        % Figure in upper right corner
        %figure('Position',[scrsz(3)/2-5 scrsz(4)/2-73 scrsz(3)/2 scrsz(4)/2]);
        % Figure fullscreen
        figure('Position',[0 0 scrsz(3) scrsz(4)]);
        
        
        subplot(2,1,1);
        hh = surf(u);
        
        set(hh,'facecolor','texturemap','edgecolor','none');
        %shading interp
        
        caxis([-0.05*maxkhsqft 0.2*maxkhsqft])
        colormap bone
        %brighten(0.3)
        %colorbar
        axis ([0,Nx+1,0,Ny+1]);
        axis equal
        axis tight
        title('hello')
    end
    if mod(t,5)==0 && animate==1
        set(hh,'zdata',u)
        %if plotoption==1
            subplot(2,1,2);
            polarplot(theta,abs(outputsflat(t,1:numouts)),'--r')
            rlim([0,0.5*maxkhsqft])
        %end
        drawnow limitrate
    end
    
end
toc(tstart)

for thetaOUT=1:numouts
    maxoutflat(thetaOUT,pass) = abs(max(outputsflat(:,thetaOUT)));
    energyflat(thetaOUT,pass) = (norm(outputsflat(:,thetaOUT).^2))/length(outputsflat(:,thetaOUT));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REDO Boundry Conditions (to take lid off diffusor)
% Matrix J: true=air, false=boundry
% Matrix K: Weighting based on J used in wave equation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% lid off diffusor
J(d+1,b1:b2) = 1;

% left and right most ends of diffuser
b1 = div(1)-round(w/2);
b2 = div(S*nW)+round(w/2);

% dividers
J(2:d+1,div(1:S*nW)) = 0;
% well bottoms
for n=2:S*nW
    J(bot(n)+1,div(n-1):div(n)) = 0;
end

% weighting matrix for wave equation
K = zeros(Ny+1,Nx+1);
K(2:Ny,2:Nx) = 2-l_sq*(J(3:Ny+1,2:Nx)+J(1:Ny-1,2:Nx)...
    +J(2:Ny,3:Nx+1)+J(2:Ny,1:Nx-1));

%update input frequency
if input=='Sine';
    f = sin(2*pi*freqIN(pass)*time);
    f((length(f)+1-winlen):length(f)) = halfHann(1:winlen).*f((length(f)+1-winlen):length(f));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% main loop - with diffusor
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
u2 = zeros(Ny+1,Nx+1);
u1 = u2;
u  = u1;

tstart = tic;
for t=1:T
    
    % 2D wave equation
    u(2:Ny,2:Nx) = l_sq *( u1(3:Ny+1,2:Nx) + u1(1:Ny-1,2:Nx) ...
        + u1(2:Ny,3:Nx+1) + u1(2:Ny,1:Nx-1))...
        + u1(2:Ny,2:Nx).*(K(2:Ny,2:Nx))     ...
        - u2(2:Ny,2:Nx);
    
    % absorb wave at bottom wall either side of diffuser
    u(2,2:div(1)) = fac*(2*(1-l_sq)*u1(2,2:div(1)) ...
        + 2*l_sq*(u1(3,2:div(1)))...
        + lambda*u2(2,2:div(1)) ...
        - u2(2,2:div(1)));
    u(2,div(S*nW):Nx) = fac*(2*(1-l_sq)*u1(2,div(S*nW):Nx) ...
        + 2*l_sq*(u1(3,div(S*nW):Nx))...
        + lambda*u2(2,div(S*nW):Nx) ...
        - u2(2,div(S*nW):Nx));
    
    % absorb wave at top wall
    u(Ny,:) = fac*(2*(1-l_sq)*u1(Ny,:) + 2*l_sq*(u1(Ny-1,:))...
        + lambda*u2(Ny,:) - u2(Ny,:));
    
    % turn on mics and absorb on sides after input bounces of diffusor
    if t > bounce
        % absorb wave at left wall
        u(:,2) = fac*(2*(1-l_sq)*u1(:,2) + 2*l_sq*(u1(:,3))...
            + lambda*u2(:,2) - u2(:,2));
        
        % absorb wave at right wall
        u(:,Nx) = fac*(2*(1-l_sq)*u1(:,Nx) + 2*l_sq*(u1(:,Nx-1))...
            + lambda*u2(:,Nx) - u2(:,Nx));
        
        
        % output
        for thetaOUT=1:numouts
            outputs(t,thetaOUT) = u(abs(round(r*sin(theta(thetaOUT))+d)),round(r*cos(theta(thetaOUT))+Q));
        end
    end
    
    % zero value at boundries
    for n=b1:b2
        for m=2:d+1
            if J(m,n)==0
                u(m,n)=0;
            end
        end
    end
    
    % input function
    if t <= length(f)
        u(P,:) = khsq*f(t);
    end
    
    % advance previous two time steps
    u2 = u1;
    u1 = u;
     
    % draw
    if mod(t,5)==0 && animate==1
        set(hh,'zdata',u)
        subplot(2,1,2);
        polarplot(theta,abs(outputs(t,1:numouts)),'--r')
        rlim([0,0.5*maxkhsqft])
        drawnow limitrate
    end
    
end
toc(tstart)

for thetaOUT=1:numouts
    maxout(thetaOUT,pass) = abs(max(outputs(:,thetaOUT)));
    energy(thetaOUT,pass) = (norm(outputs(:,thetaOUT).^2))/length(outputs(:,thetaOUT));
end

bigger(pass) = max(maxoutflat(:,pass));
if bigger(pass)<max(maxout(:,pass))
    bigger(pass)=max(maxout(:,pass));
end

'countdown'
repeat-pass

end

big=zeros(repeat,1);
bigflat=zeros(repeat,1);
for pass=1:repeat
    big(pass) = max(maxout(:,pass));
    bigflat(pass) = max(maxoutflat(:,pass));
end
biggest = max(big);
if biggest<max(bigflat)
    biggest=max(bigflat);
end

strnW = num2str(nW);
strS = num2str(S);
strlo = num2str(f_lo);
strhi = num2str(f_hi);
titlee = strcat('Bandwidth: ',strlo,'Hz-',strhi,'Hz','_nW=',strnW,'_S=',strS);
scrsz = get(0,'ScreenSize');
fig = figure('Position',[0 0 scrsz(3) scrsz(4)]);     
fig.Name= titlee;
for pass=1:repeat
subplot(ceil(repeat/3),ceil(repeat/3),pass);
polarplot(theta,-maxoutflat(:,pass),'--k',theta,-maxout(:,pass),'-r')
strfreq = num2str(freqIN(pass));
strFREQ = strcat(strfreq,'Hz');
title(strFREQ),rlim([0,biggest])
end

