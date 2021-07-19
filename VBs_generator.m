%% Vector beams obtained as a combination of polarized high-order Gaussian modes
% Coded by David Marco Castillo, Universidad Miguel Hernández de Elche,
% dmarco@umh.es
% Copyright (c) 2021, David Marco Castillo.
% Last update: 18/07/2021

% This code may be freely used, modified and
% distributed under the GNU General Public License: 
% http://www.gnu.org/copyleft/gpl.html
%Not for commercial use.

%The code uses the functions HermitePoly and LaguerreL coded by David 
%Terr and Manuel Diaz, respectively. These functions calculate the Hermite 
%and generalized Laguerre polynomials required for the Hermite-Gaussian and
%Laguerre-Gaussian beams. The two functions are attached as separated 
%files. I do not own the rights of these two functions.

%************************************************************************%


%This code plots the intensity profile and the polarization map of a vector
%beam obtained as the superposition of two polarized high-order Gaussian modes. It
%also plots the intensity profile and the phase of the combined modes. The
%software can be extended to arbitrary combinations of modes at each
%polarization component by creating new combinations in the function
%VectorBeam.

%The generated vector beam can be changed by varying parameters in sections
%1 to 4, including the added modes, their polarization, their amplitude and
%phase relations and the plot parameters. For deeper changes you will have 
%to go inside the functions defined below (after section 7) or to the 
%section 7, where almost all the code is executed. In addition, an 
%elliptical polarizer (section 5) and and/or an elliptical retarder
%(section 6) can be placed after the vector beam.

close all
clear
%% 1. Beam parameters
% General Gaussian beam parameters. All the Gaussian modes share the same
% parameters. You can make each Gaussian beam to have different parameters
% by changing the z, wavelength and w0 parameters that you send to each
% LaguerreGauss or HermiteGauss functions in the function VectorBeam.
z = 0; % z coordinate in meters
wavelength = 633e-9; % wavelength in meters
w0 = 0.0015; % waist in meters
%% 2. Polarization vectors
% orthonormal vectors in the linear |x> |y> basis
polcomponentA = (1/sqrt(2)) * [1 1i]; %1st component for mode A
polcomponentB = (1/sqrt(2)) * [1 -1i]; %2nd component for mode B
%right and left circular polarizations by default.
%% 3. Type of mode addition
%Choose the type of mode at each vector component:
modesumtype =1;
    %1. modesumtype =1 -> Laguerre-Gaussian modes sum.
    %2. modesumtype =2 -> Hermite-Gaussian modes sum.
    %3. modesumtype =3 -> Custom mode sum.   
%1: Laguerre-Gaussian (LG) modes sum. One LG mode at each vector component
%mode A:
pA = 0; %radial number
lA = -3; %azimuthal number
%mode B:
pB = 2; %radial number
lB = 1; %azimuthal number

%2: Hermite-Gaussian (HG) modes sum. One HG mode at each vector component.
%modeA:
mA = 1; 
nA = 0;
%modeB:
mB = 0; 
nB = 1;

%For modesumtype = 1, 2:
% Phase and amplitude relations between modes at each polarization
% component:
%VectorBeam = 
%cos(modeamplitudefactor)*exp(-1i*modephasefactor)*modeA*|polcomponentA> +
%sin(modeamplitudefactor)* exp(1i*modephasefactor)*modeB*|polcomponentB>
polcomponentsamplitudefactor = pi/4; % Amplitude relation between polarization components
polcomponentsphasefactor = 0; % Phase relation between polarization components


%3: Custom sum. Go down to the function VectorBeam and customize the mode
%addition. Mix LG and HG modes at each vector component with the amplitude
% and phase relations at will.


%% 4. Plot parameters

rangescalingfactor = 1.5; %defines the plot range for the transverse section of the beam.
%For rangescalingfactor = 1, the plot range is the diameter of the beam 2w,
%where w is the width of the fundamental Gaussian mode:
%w=w0*sqrt(1+(z/zR)^2).


%%4.1. Intensity and phase parameters 

showmodeproperties = true; % set true to plot the intensity and phase of the
%modes at each polarization component
plotintensity = true; %plot intensity behind the vector beam polarization map
%Set to 0 to only see the polarization map, set to 1 to see the intensity and
%the polarization map superimposed.
% Resolution for the Intensity and phase
resolutionintensityandphase =1000; % number of points per Row/Column for
%plotting the intensity of each mode and of the resulting vector beam and 
%the phase of each mode.

%Intensity interpolarion:
interpolateintensity = true; % set true to interpolate the intensity plots
numberofpointsinterpolatedintensity = 2000;
intensitytransparency = 1; %set the transparency of the intensity plot.



%%4.2. Polarization Map parameters

plotpolarizationmap = true; % set to true to plot the polarization map
sizeEllipsefactor = 0.031; % 0.031 gives a nice depiction of the map for 19 ellipses per row/colum.
NEllipses = 19; % number of ellipses per Row/Column.
linearthreshold = 0.02; % threshold for phase and amplitude to consider
%a polarization state as linear.
intensitythresholddrawellipse = 0.01; % if the normalized intensity
%of the beam is lower than this value no ellipse will be plotted at that
%point.

%Ellipses visual properties:
ellipselinethickness = 1.6;
linearellipsecolor = [0 1 0]; % Linear polarization: Green
righthandedellipsecolor = [0 0 1]; % Right-handed polarization: Blue
lefthandedellipsecolor = [1 0 0]; % Left-handed polarization: Red
%% 5. Place an elliptical retarder after the vector beam
switchonWaveplate = false; % switchonWaveplate = true places an elliptical 
%retarder after the vector beam
if switchonWaveplate == true
    eigenaxesorientation = 0;
    retardance = pi/2;
end
%% 6. Place an elliptical polarizer after the vector beam
switchonPolarizer = false; % switchonPolarizer = true places an elliptical 
%polarizer after the vector beam.
%Polarizer parameters 
if switchonPolarizer == true
    eigenstateorientation =pi/4;
    eigenstatephasecomponents =pi/2; % phase between the  
% components of the eigenstate vector that defines the polarizer.
end

%If both waveplate and polarizer are switched on, the waveplate will act first.
%-------------------------------------------------------------------------------------------------------%
%End of the customizable parameters section
%% 7. Plot Section %%
zR=pi*w0^2/wavelength; %Rayleight length
w=w0*sqrt(1+(z/zR)^2); %Radius of the beam
range = (w*2)*rangescalingfactor; % range selected for X and Y variables (diameter 
%of the beam * scaling factor).

if  plotintensity == true || showmodeproperties == true
    [x,y,X,Y,R,THETA] = Generatepointsinspace(range,resolutionintensityandphase);   
    [modeA,modeB,Jx,Jy] = VectorBeam(X,Y,R,THETA,z,wavelength,w0,pA,pB,lA,lB,mA,mB,nA,nB,modesumtype,polcomponentsamplitudefactor,polcomponentsphasefactor,polcomponentA,polcomponentB);
    %Plot amplitude and phase of the modes
    if showmodeproperties == true
        figure,
        PlotModeIntensityandPhase (range,x,y,X,Y,modeA,modeB,interpolateintensity,numberofpointsinterpolatedintensity)
    end

figure, 
    %Plot vector beam intensity
    if plotintensity == true
        PlotIntensity(range,x,y,X,Y,Jx,Jy,interpolateintensity,numberofpointsinterpolatedintensity)
        colormap bone
        colorbar
        alpha(intensitytransparency)
        datacursormode on
    end
end

%Plot polarization map
if plotpolarizationmap == true
    sizeEllipse = range*sizeEllipsefactor; % size of ellipses
    [xPolMap,yPolMap,XPolMap,YPolMap,RPolMap,THETAPolMap] = Generatepointsinspace(range,NEllipses);
    [modeA,modeB,JxPolMap,JyPolMap] = VectorBeam(XPolMap,YPolMap,RPolMap,THETAPolMap,z,wavelength,w0,pA,pB,lA,lB,mA,mB,nA,nB,modesumtype,polcomponentsamplitudefactor,polcomponentsphasefactor,polcomponentA,polcomponentB);
    PlotPolarizationMap(JxPolMap,JyPolMap,NEllipses,range,sizeEllipse,linearthreshold,intensitythresholddrawellipse,ellipselinethickness,linearellipsecolor,righthandedellipsecolor,lefthandedellipsecolor)
end

%% Waveplate action
if switchonWaveplate == true
    [Jx,Jy] = Waveplate(eigenaxesorientation, retardance,Jx,Jy);
    [JxPolMap,JyPolMap] = Waveplate(eigenaxesorientation, retardance,JxPolMap,JyPolMap);
end
%% Polarizer action
if switchonPolarizer == true
    [Jx,Jy] = Polarizer(eigenstateorientation, eigenstatephasecomponents,Jx,Jy);
    [JxPolMap,JyPolMap] = Polarizer(eigenstateorientation, eigenstatephasecomponents,JxPolMap,JyPolMap);
end

%% Functions

function [x,y,X,Y,R,THETA] = Generatepointsinspace(range,resolution)
    x = linspace(-range, range,resolution);
    y = linspace(range,-range,resolution);
    [X,Y] = meshgrid(x,y);
    
    [THETA,R] = cart2pol(X,Y);
    THETA = THETA - 2*pi*floor(THETA/(2*pi));
end

function PlotIntensity(range,x,y,X,Y,Jx,Jy,interpolate,numberofinterpolationpoints)
    intensity = abs(Jx).^2 + abs(Jy).^2;
    intensitynormalized = intensity/(max(intensity(:)));



    if interpolate == true
        xq = linspace(-range, range,numberofinterpolationpoints);
        yq = linspace(range, -range,numberofinterpolationpoints);
        [Xq,Yq] = meshgrid(xq,yq);
        intensitynormalizedfinal = interp2(X,Y,intensitynormalized,Xq,Yq,'spline');
    else
        xq=x;
        yq=y;
        intensitynormalizedfinal = intensitynormalized;
    end
    imagesc(xq,yq,intensitynormalizedfinal)
    set(gca,'YDir','normal')
    xlabel('x');
    ylabel('y');
    caxis([0 1])
    pbaspect([1 1 1])
end

function PlotPolarizationMap(Jx,Jy,NEllipses,range,sizeEllipse,tol,intensitythresholddrawellipse,ellipselinethickness,linearellipsecolor,righthandedellipsecolor,lefthandedellipsecolor)
    jump = 2*range/(NEllipses-1); %distance between ellipses
    modJout = sqrt((abs(Jx).^2)+(abs(Jy).^2));
 
    modJout2 = modJout.^2; %intensity of the ellipses
    Ex0 = abs(Jx)./modJout;
    Ey0 = abs(Jy)./modJout;
    %Phases
    deltax = angle(Jx);
    deltay = angle(Jy);
    Delta = deltay-deltax; %Overall phase for the y component.  
    t=0:pi/500:2*pi; %parameter for plotting the ellipse
    
    hold on
    for i = 1:NEllipses
        for j = 1:NEllipses
            if modJout2(i,j)/max(modJout2(:)) > intensitythresholddrawellipse %%if the intensity
            %is large enough, plot the ellipse.
                    
            %Draw the ellipses in their corresponding position in
            %the map:
            Ex = Ex0(i,j)*sin(t)*sizeEllipse - range + (j-1)*jump;
            Ey = Ey0(i,j)*sin(t+Delta(i,j))*sizeEllipse + range - (i-1)*jump;

            
            %%Ellipse Handeness:           
            d = Delta(i,j);            
            %Linear polarization
            if abs(d) < tol || abs(abs(d) - pi) < tol || abs(abs(d) - 2*pi) < tol || abs(Ex0(i,j)) < tol || abs(Ey0(i,j)) < tol
                plot(Ex,Ey,'color',linearellipsecolor,'LineWidth',ellipselinethickness)
            %Right-handed polarization
            elseif d > 0 && d <pi || d > -2*pi && d <-pi
                plot(Ex,Ey,'color',righthandedellipsecolor,'LineWidth',ellipselinethickness)
            %Left-handed polarization
            else
                plot(Ex,Ey,'color',lefthandedellipsecolor,'LineWidth',ellipselinethickness)
            end
            pbaspect([1 1 1])
          
            else
            end
       
        end %end for loop draw ellipses at rows
    end %end for loop draw ellipses at columns

end

function [modeA,modeB,Jx,Jy] = VectorBeam(X,Y,R,THETA,z,wavelength,w0,pA,pB,lA,lB,mA,mB,nA,nB,modetype,modeamplitudefactor,modephasefactor,e1,e2)
    % Addition of LG modes
    if modetype == 1
        modeA =  LaguerreGauss(z,wavelength,w0,pA,lA,R,THETA);
        modeB =  LaguerreGauss(z,wavelength,w0,pB,lB,R,THETA);
    end

    %%Addition of HG modes
    if modetype == 2
        modeA = HermiteGauss (z,wavelength,w0,mA,nA,X,Y,R);
        modeB =  HermiteGauss (z,wavelength,w0,mB,nB,X,Y,R);
    end
    
    %%Free addition of LG and HG modes
    if modetype == 3
    %Feel free to add the modes you want at each vector component with
    %any phase and amplitude relations:
        modeA = LaguerreGauss(z, wavelength,w0,0,-4, R, THETA)+ LaguerreGauss(z, wavelength,w0,0,8, R, THETA);
        modeB = LaguerreGauss(z,wavelength,w0,0,-1,R,THETA);
    end


    %Modes at each vector component
    modee1 = cos(modeamplitudefactor)*exp(-1i*modephasefactor/2)*modeA;
    modee2 = sin(modeamplitudefactor)* exp(1i*modephasefactor/2)*modeB;

    % Output Jones Vector
    Jx = modee1 * e1(1) + modee2 * e2(1);
    Jy = modee1 * e1(2) + modee2 * e2(2);
end

function PlotModeIntensityandPhase (range,x,y,X,Y,modeA,modeB,interpolateintensity,numberofpointsinterpolatedintensity)

        %%Mode A
        %Intensity Mode A
        modeproperties(1) = subplot(2,2,1);
        
        JxA = abs(modeA);
        JyA = 0;
        PlotIntensity(range,x,y,X,Y,JxA,JyA,interpolateintensity,numberofpointsinterpolatedintensity)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        colormap(modeproperties(1),hot)
        colorbar
        title('mode A intensity')

        %Phase Mode A
        modeproperties(2) = subplot(2,2,2);
        modeAphase = angle(modeA);
        modeAphase = modeAphase - 2*pi*floor(modeAphase/(2*pi));

        imagesc(x,y,modeAphase);
        caxis([0 2*pi])
        colormap(modeproperties(2),jet)
        colorbar
        pbaspect([1 1 1])
        set(gca,'YDir','normal')
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        title('mode A phase')
        xlabel('x');
        ylabel('y');
        
        %%Mode B
        %Intensity Mode B        
        modeproperties(3) = subplot(2,2,3);
        JxB = abs(modeB);
        JyB = 0;
        PlotIntensity(range,x,y,X,Y,JxB,JyB,interpolateintensity,numberofpointsinterpolatedintensity)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        colormap(modeproperties(3),hot)
        colorbar
        title('mode B intensity')
        
        %Phase Mode B
        modeproperties(4) = subplot(2,2,4);
        modeBphase = angle(modeB);
        
        modeBphase = modeBphase - 2*pi*floor(modeBphase/(2*pi));

        imagesc(x,y,modeBphase);
        caxis([0 2*pi])
        colormap(modeproperties(4),jet)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        set(gca,'YDir','normal')
        colorbar
        pbaspect([1 1 1])
        title('mode B phase')
        xlabel('x');
        ylabel('y');
        hold off
end

function [LaguerreGaussBeam] = LaguerreGauss(z,wavelength,w0,p,l,R,THETA)
%Returns a complex 2-dimensional matrix. Each element of the matrix contains
%the amplitude and phase of a Laguerre Gaussian beam at a coordinate (x,y) 
%in space given as a complex number.

%Parameters:
    %*z: z coordinate in space
    %*wavelength
    %*w0: waist
    %*p: radial index
    %*l: azimuthal index
    %*R: 2-dimensional matrix with the info about the radial coordinate for each element
    %*THETA: 2-dimensional matrix with the info about the polar coordinate for each element

A=factorial(p)*sqrt(2/(pi*factorial(p)*factorial(abs(l)+p))); %Amplitude term normalized
zR=pi*w0^2/wavelength; %Rayleight length
Rc=z+zR^2/z; %Radius of curvature
w=w0*sqrt(1+(z/zR)^2); %Radius of the beam
k=2*pi/wavelength; %wave number
gouyphase=((2*p+abs(l))+1)*atan(z/zR); %Gouy phase


Amplitude=(A/w)*((sqrt(2)*R/w).^(abs(l))).*LaguerreL(p,abs(l),2*(R.^2)/w^2).*exp(-R.^2/w^2);
Phase = (k*(R.^2)/(2*Rc))+l*THETA-gouyphase;
Phase= Phase - 2*pi*floor(Phase/(2*pi));
LaguerreGaussBeam = Amplitude.*exp(1i*Phase);
end

function [HermiteGaussBeam] = HermiteGauss (z,wavelength,w0,m,n,X,Y,R)
%Returns a complex 2-dimensional matrix. Each element of the matrix contains
%the amplitude and phase of a Hermite Gaussian beam at a coordinate (x,y) 
%in space given as a complex number.

%Parameters:
    %z: z coordinate in space
    %wavelength
    %w0: waist
    %m and n: mode indices
    %X and Y: 2-dimensional matrices with the info about the x and y coordinates for each element
    %R: 2-dimensional matrix with the info about the radial coordinate for each element
    
k=2*pi/wavelength;
A=sqrt((2^(1-(n+m)))/(pi*factorial(n)*factorial(m)));
zR=pi*w0^2/wavelength;
Rc=z+zR^2/z;
gouyphase=((n+m)+1)*atan(z/zR);
w=w0*sqrt(1+(z/zR)^2);

Amplitude=(A/w)*polyval(HermitePoly(m),sqrt(2)* X/w).*polyval(HermitePoly(n),sqrt(2)* Y/w).*exp(-(R.^2)/w^2);
phase= exp(1i*k*(R.^2)/(2*Rc)).*exp(-1i*gouyphase);

HermiteGaussBeam = Amplitude.*phase;

end

function [Joutx,Jouty] = Waveplate(eigenaxesorientation, retardance,Jinx,Jiny)
    %Returns the ouput polarization state after passing throught a
    %waveplate.
    %Returns 2 complex matrices. One with the |x> component of the Jones vector at each point
    %in space and other with the |y> component.
    
    %Parameters:
    %*eigenaxesorientation: orientation of the eigenaxes of
    %the waveplate in radians.
    %*retardance: retardance of the waveplate in radians.
    %*Jinx: |x> component of the input Jones vector at each point in space
    %*Jiny: |y> component of the input Jones vector at each point in space
    
    Joutx = (cos(eigenaxesorientation)^2 * exp(1i*retardance/2) + sin(eigenaxesorientation)^2 * exp(-1i*retardance/2)) * Jinx + cos(eigenaxesorientation) * sin(eigenaxesorientation) * 2 * 1i * sin(retardance/2) * Jiny;
    Jouty = cos(eigenaxesorientation) * sin(eigenaxesorientation) * 2 * 1i * sin(retardance/2) * Jinx + (sin(eigenaxesorientation)^2 * exp(1i*retardance/2) + cos(eigenaxesorientation)^2 * exp(-1i*retardance/2)) * Jiny;
end

function [Joutx,Jouty] = Polarizer(eigenstateorientation, eigenstatephasecomponents,Jx,Jy)
    %Returns the ouput polarization state after passing throught a general
    %elliptical polarizer.
    %Returns 2 complex matrices. One with the |x> component of the Jones vector at each point
    %in space and other with the |y> component.
    
    %Parameters:
    %*eigenstateorientation: orientation of the transmission axis of the
    %polarizer in radians.
    %*eigenstatephasecomponents: %Phase between the components of the 
    %eigenstate vector that defines the polarizer
    %*Jinx: |x> component of the input Jones vector at each point in space
    %*Jiny: |y> component of the input Jones vector at each point in space
    
    
    Joutx = cos(eigenstateorientation)^2 * Jx + exp(-1i*eigenstatephasecomponents) * sin(eigenstateorientation) * cos(eigenstateorientation) * Jy;
    Jouty = Jx * exp(1i*eigenstatephasecomponents) * sin(eigenstateorientation) * cos(eigenstateorientation) + sin(eigenstateorientation)^2 * Jy;
end