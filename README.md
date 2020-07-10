# OCT_MC  
Monte-Carlo simulation of an FD-OCT B-Scan.  
  
Description  
  
===Missing Features / Bugs (to be implemented)===  
-Check why a low value of mua and g (maybe mus) slows the simulation a lot  
-Need to change how the BScan width is controlled  
-Making a map of the code  
-Testing the functions  
-Removing the extra code  
-Fresnel reflection  
-Mirror reflection  
-Refraction  
-Check if Noise is good  
-Gaussian beam source  
-Wavelength width on the detection  
-Focusing/diverging beams instead of perfectly foccus  
  
===Added Features===  
-Electric filter on the OCT raw A-line 
-The thesis proposed mesh intersection is the same as in St-Jacques code.  
-Corrected the bug in the difference in size from output and input. It was due to the input file not saving all decimals on the dx/dy/dz parameters. 0.00045 became 0.0004.  
-Selecting a_coef as a parameter  
-Fix progress bar of the .c simulation  
-Selecting the number of simulating photons as oppose to time-min.  
-Selecting the number of A-line in the matlab script.  
-Selecting the parameter p  
-Photon path length are now counted in each individual medium.  
-OCT matlab script now allows for taking into account the refractive index impact on the OCT image.  
-Flat noise was added  
  
MergedCode.c  
  
===Parameters===  
a_coef: Biasing coefficient of the importance sampling. Default = 0.5  
p: Probability of additional bias scattering. Default = 0.925  
  
  
===Functions===  
RandomGen: Generates a double precision random number ranging from 0 to 1 inclusive. Based on W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, "Numerical Recipes in C," Cambridge University Press, 2nd edition, (1992).  
RandomNum: Shorten version of RandomGen to be more convient. RandomNum = (double) RandomGen(1, 0, NULL)  
SameVoxel: Determine if the two position are located in the same voxel. Returns 1 if same voxel, 0 if not same voxel.  
max2: Fast calculation of which between 2 values is the biggest. (checked)  
min2: Fast calculation of which between 2 values is the smallest. (checked)  
min3: Fast calculation of which between 3 values is the smallest. (checked)  
RFresnel: Computes reflectance as photon passes from medium 1 to medium 2 with refractive indices n1,n2. Incident angle a1 is specified by cosine value ca1 = cos(a1). Program returns value of transmitted angle a1 as value in *ca2_Ptr = cos(a2).  
FindVoxelFace: (not used) the boundary is the face of some voxel. Find the the photon's hitting position on the nearest face of the voxel and update the step size.
FindVoxelFace2: Modified version of Jacques code.  
  
===Variables===  
a: Variable used in the original code to determine if two positions are in the same voxel  
b: Variable used in the original code to determine if two positions are in the same voxel  
c: Variable used in the original code to determine if two positions are in the same voxel  
W: Photon weigth (max(W) = 1)  
photon_status: ALIVE or DEAD? Might be only partially dead if a split occured. (The photon died but his split survived.)  
i_photon: increment of photon photon packet  
c_photon: count of collected photon packet  
j: increment of a for loop. used twice?  
NN: Number of voxel in the 3D simulated volume  
Nyx: Number of voxel in a 2D cut on the y and x axis of the simulated volume  
Nx,Ny,Nz: Size of the simualted space in voxel in the corresponding direction  
dx,dy,dz: Width of a voxel in cm  
xs,ys,zs: Control of the position of the source in cm? xs and ys = 0 for a centered source.  
x,y,z: Current position of the photon  
x_cont,y_cont,z_cont: x,y,z coordonates of the continuing photon that splited  
tempx,tempy,tempz: Temporary value for x,y,z.  
ux,uy,uz: Direction of the photon before the scattering step  
uxx, uyy, uzz: Direction of the photon just after scattering for the continuing photon that splited.  
ux_cont,uy_cont,uz_cont: ux,uy,uz coordonates of the continuing photon that splited  
ix,iy,iz: Indexes used to determine if the photon goes through a border and change medium  
vx,vy,vz: Unit vector pointing from the current position of the photon toward the direction of the selected detector  
upx, upy, upz: Direction of the photon just after scattering, biased.  
detx: Position of the dectector/launch position. Ranges from -radius to +radius  
dety: Position of the detector on the y axis. Centered with dety = 0.  
detz: Position of the detector on the z axis. detz = det_z  
det_z: Position of the dectector in z. Redundant, always = zs.  
F: Fluence rate  
R: Reflectance It was specify in the original code that it was not ready to be used. We don't use it either.  
Rd: I think it quantifies the total quantity of photon reflected at the surface. Not used.  
det_num: Detector ID. Can have a value of -1 if no detection has occured yet.  
first_bias_done: 0 or 1? Determine if the photon has been bias or not yet.  
cont_exist: 0 or 1? determine if the photon have been split yet.  
L_current: Photon likelihood. Initially = 1  
L_cont: L_current of the continuing photon after a split  
s_total: photon total path length. Starts at 0.  
s_total_cont: s_total of the continuing photon that splited?  
z_max: photon depth reached  
z_max_cont: z_max of the continuing photon that splited?  
Ndetectors: Number of detectors/Alines  
det_radius: Radius of the detector in cm  
cos_accept: Aperture of the detector/cos of the angle  
Pick_det: Picked detector. Ranged from 1 to Ndetectors in   
radius: NOT A RADIUS. It is half the B-line width.  
detx: Position in cm? of the current A-line position  
launchflag: if = 1 Manually set launch direction to straigth down.  
mcflag: 0 = uniform beam. Other values are for 1 = gaussian, 2 = isotropic and 3 rectangular beam  
bflag: 1 = The photon is still inside the simulated space. 0 = it no longer is.  
surfflag: 0 = photon is inside the tissue part of the simulated space. 1 = photon is inside the air part of the simulated space.  
boundaryflag: Determines where the photon can escape. 0 = no boundaries, 1 = escape at all boundaries, 2 = escape at surface only From the original code.  
a: Random variable ranging from 0 to 1 and including 0 and 1  
rnd: Random variable ranging from 0 to 1 sometimes including or excluding one or both of the extremities  
sleft: Distance the photon have left to propagate but without considering mus. (Better formulate)  
s: Random propagating distance in cm  
sv: 1 = The photon is in the same voxel as before. 0 = The photon have changed voxel.  
absorb: Proportion of the photon packet absorb  
W: Weigth of the photon  
W_cont: Weigth of the photon of the continuing photon that splited?  
i: Index of the voxel. All voxels in the simulated space have a specific index associated to it.  
i_cont: i of the continuing photon that splited? The contuing photon is the one that didn<t scattered.  
ls: A tiny step value to ensure that the photon is not exactly between two voxels. In cm.  
v: Related to the tissue type  
type: Current tissue type where the photon is located  
mua: Current absorption coeficient in cm  
muav: Possible absorption coeficients in the simulation in cm  
mus: Current scattering coeficient in cm  
musv: Possible scattering coeficients in the simulation in cm  
g: Current anisotropy  
gv: Possible anisotropy in the simulation  
cont_exist: Determines if a photon will split and continue moving? 0=no and 1=yes  
first_bias_done: Determines if a scaterred bias has occured yet =1 or not =0  
temp: Temporary variable used to calculate part of an equation in order tos shorten it.  
costheta_B: Random result of cos(theta_B). theta_B is the random bias axial angle between the detector and scattered direction  
sintheta_B: Random result of sin(theta_B). theta_B is the random bias axial angle between the detector and scattered direction  
psi: Random psi angle determining the lateral direction  
cospsi: cos(psi)  
sinpsi: sin(psi)  
costheta_S: Random result of cos(theta_B). theta_S is the random bias axial angle between the previous photon direction and scattered direction  
L_temp: 1 minus the Likelihood ratio of the current scattering event  
ONE_MINUS_COSZERO: One very small value. I think it used to avoid bugs when a value is too close to 0.  
a_coef: Biasing coefficient of the importance sampling  
p: Probability of additional bias scattering.  
f_HG: Result of the Henyey Greenstein function  
f_B: Result of the bias function  
CHANCE: Chance of the photon dying when under the weigth. Defined as 0.1  
Ntiss: Maximum number of tissue types. (Can be edited to be bigger)  
DetS: Path Length of the detected photons  
DetW: Weigth of the detected photons  
DetL: Likelihood ratios of detected photons  
DetID: IDs of detected photons (range from 1 to 512 for 512 detectors/emittors)  
c_photon: Current photon packet itteration?  
m: Indices used in for loops.  
