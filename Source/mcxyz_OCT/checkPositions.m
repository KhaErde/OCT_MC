myname = '6_1_Compare';

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%Number of photon loaded at the same time
nph = 200000; %%%%CHOOSE VALUE%%%%
%Limit the total number of loaded photons if you don't want to analyse all the data
%Input 0 if you want to load the whole data
maxnph = 0;
%OCT wavelengths IN CENTIMETERS
lambda_start = 1150e-7; %%%%CHOOSE VALUE%%%%
lambda_stop = 1450e-7; %%%%CHOOSE VALUE%%%%
%Number of sample point of the OCT wavelength width
samplePoints= 2048; %%%%CHOOSE VALUE%%%%
%Choose to apply electric filter. Put a high value if none
maxDepth = 0; %%%%CHOOSE VALUE%%%% 0 if none
%Compress image for refractive index
n_cor = 1;
%Chosse refractive of the different mediums. Index. Can be a function of the wavelength
rn = ones(samplePoints,1); %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,1) = 1; %%%%CHOOSE VALUE%%%%
%rn(1:samplePoints,2) = 1.4; %%%%CHOOSE VALUE%%%%
%rn(1:samplePoints,3) = 1; %%%%CHOOSE VALUE%%%%
%rn(1:samplePoints,4) = 1; %%%%CHOOSE VALUE%%%%

%Choose amplitude of noise
noise_amp = 0; %%%%CHOOSE VALUE%%%%

%% parameters
n = 1;
time_min = A(n); n = n + 1;
a_coef = A(n); n = n + 1;
p = A(n); n = n + 1;
Ndetectors = A(n); n = n + 1;
det_radius = A(n); n = n + 1;
cos_accept = A(n); n = n + 1;
lambda = A(n); n = n + 1;
f = A(n); n = n + 1;
D = A(n); n = n + 1;
z_f_img = A(n); n = n + 1;
h_step = A(n); n = n + 1;
Nx = A(n); n = n + 1;
Ny = A(n); n = n + 1;
Nz = A(n); n = n + 1;
dx = A(n); n = n + 1;
dy = A(n); n = n + 1;
dz = A(n); n = n + 1;
mcflag = A(n); n = n + 1;
launchflag = A(n); n = n + 1;
boundaryflag = A(n); n = n + 1;
xs = A(n); n = n + 1;
ys = A(n); n = n + 1;
zs = A(n); n = n + 1;
xfocus = A(n); n = n + 1;
yfocus = A(n); n = n + 1;
zfocus = A(n); n = n + 1;
ux0 = A(n); n = n + 1;
uy0 = A(n); n = n + 1;
uz0 = A(n); n = n + 1;
radius = A(n); n = n + 1;
waist = A(n); n = n + 1;
zsurf = A(n); n = n + 1;
x_start = A(n); n = n + 1;
y_start = A(n); n = n + 1;
z_start = A(n); n = n + 1;
Nt = A(n);  n = n + 1;
j = n;
for i=1:Nt
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
    j=j+1;
    nrv(i,1) = A(j);
    j=j+1;
end

rn(1:samplePoints,1) = nrv; %%%%CHOOSE VALUE%%%%

%% Load x position of photon
filename = sprintf('%s_Track_x.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    Track_x = fread(fid, 'double');
    fclose(fid);
toc

%% Load y position of photon
filename = sprintf('%s_Track_y.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    Track_y = fread(fid, 'double');
    fclose(fid);
toc
%% Load z position of photon
filename = sprintf('%s_Track_z.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    Track_z = fread(fid, 'double');
    fclose(fid);
toc
%% 
%{
figure()
plot(Track_z, Track_x)
xlabel('z [cm]')
ylabel('x in [cm]')
title("z vs x")

figure()
plot(Track_z, Track_y)
xlabel('z [cm]')
ylabel('y in [cm]')
title("z vs y")
%}
figure()
Track_xy = sqrt(Track_x.^2 + Track_y.^2);
plot(Track_z, Track_xy,'LineStyle',':')
xlabel('z [cm]')
ylabel('sqrt(x^2+y^2) in [cm]')
title("z vs sqrt(x^2+y^2)")
txt = ['Focal length = ' num2str(f) 'cm, Beam diameter = ' num2str(D) 'cm, imaging lens position = ' num2str(z_f_img) 'cm, n = ' num2str(nrv)];
subtitle(txt)
hold on
%%

%Calculate beam width at the sample's surface
z_f = z_f_img + f; %Location of the focus beam in the medium for nr = 1
w_0 = (2*lambda*f)/(pi*D); %Beam radius at minimum waist
z_R = pi*w_0^2/lambda; %Rayleigth range
w_surf = w_0*sqrt(1+(z_f/z_R)^2); %Waist at the surface
%sigma_surf = w_surf/2; %Sigma of the gaussian distribution of the light 
                       %intensity at the surface
sigma_surf = 2 * w_surf;
%Generate a random number from a gaussian distribution
r = normrnd(0,sigma_surf); %You'll need to create the function in C but 
                           %multiple examples are available. Just test it
                           %before putting it in the code. Run it and look
                           %at the distribution it gives you.
theta = rand(1)*2*pi;



% ... other part of the code until the photon is actualy launch
detx = 0; %for the sake of this example



%Photon packet starting position and distance from the focal point
x = sin(theta)*r + xs + detx;
y = cos(theta)*r + ys;
z = zs;
x = x_start;
y = y_start;
z = z_start;
z_f_left = z - z_f; %Distance from the focal point if medium was vacuum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mus = 1;
nr = 1;
z_f_left_begin = z_f_left;
h_step_begin = h_step;
%Start the loop of the first photon launch before scattering
rnd = rand(1);
%sleft	= -log(rnd)*1e-5;
%sleft   = 50e-5;
sleft   = f;
%sleft = -log(rnd);
clear rnote znote
n = 1;
while sleft > 0
    %If what is left is smaller than a h_step, make the h_step smaller so
    %that the photon stops exactly at sleft.
    if h_step > sleft
        h_step = sleft;
    end
    s     = h_step; %h_step is used instead of s and mus until the first 
                    %scattering event. We are no longuer moving the photon
                    %until next scatter, but untill next ux,uy,uz
                    %approximation instead. We are still counting down from
                    %sleft to know when it scatters again.
    
    %before each time you are about to move your photon, you need to update
    %the direction
    z_R = nr*pi*w_0.^2/lambda; %Rayleigth range. Needs to be updated 
                                %incase the refractive index changed
   % R = (z_f_left.*nr).*(1+(z_R./(z_f_left.*nr)).^2); %Radius of curvature
    R = -(z_f_left*nr)*(1+(z_R/(z_f_left*nr))^2); %Radius of curvature
    % R = ((z-z_f).*nr).*(1+(z_R./(z-z_f)).^2); debug line
    temp = sqrt(1+(x^2+y^2)/R^2);         %Temporary variable to 
                                              %shorten the next 
                                              %calculations
    ux = (-x/R)/temp; %Same as in Hork 2015
    uy = (-y/R)/temp; %Same as in Hork 2015
    uz = 1/temp;
    
    rnote(n) = sqrt(x.^2+y.^2);
    znote(n) = z;
    
    %You'll need to be tempx or x as needed in the code
    
    x = x + s*ux;
    y = y + s*uy;
    z = z + s*uz;
    
    
    z_f_left = z_f_left + s*uz/nr; %Distance from the focal point if medium  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                   %was vacuum
    %z_f_left = z_f_left - s*uz/nr;   
    %... update tissue properties and everything.
    
    sleft = sleft - s*mus;
    n = n + 1;
end
plot(znote,rnote,'LineStyle','--')
legend('C-code','Matlab')
hold off
