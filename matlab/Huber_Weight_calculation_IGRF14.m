function [m_i_huber,m_median,TABLE] = Huber_Weight_calculation_IGRF14(TYPE)

% [m_i_huber,m_median,TABLE] = Huber_Weight_calculation_IGRF14(TYPE)
%
% Estimate the Robust Huber Gauss coeffcients from a set of IGRF candidate models 
% This example applies to the IGRF-14 candidate submitted October 2024 
%
% Input:
% TYPE = 'SV' for Secular Variation Models
% TYPE = 'DG' for DGRF
% TYPE = 'IG' for IGRF
%
% Output:
% m_i_huber: estimated coefficients using the Huber weight method in space.
% m_median : coefficients using the median of the candidates 
% TABLE    : A formatted table of statistical differences between Huber,
%		median and candidate models
% The code also writes a .cof file to the HOME path
%
% Requirements: 
%	Uses the Matlab Statistical Toolbox robustfit.m
%
% This approach is described and used in:
% Thébault, E., Finlay, C.C., Alken, P. et al. 
% Evaluation of candidate geomagnetic field models for IGRF-12. 
% Earth Planet Sp 67, 112 (2015) doi:10.1186/s40623-015-0273-4
% Original code from E. Thebault (October 2014)
% Last modified by C. Beggan (Nov-2024)
%
% To use this code, place the candidate model .cof files in a directory on the path, manually set (below).
% Alter the labels in LLf manually too.
% Run separately for the DGRF, IGRF and SV candidates. e.g.
% From Matlab GUI or command line, for the SV Huber weighted model:
% >>  [m_i_huber,m_median,TABLE] = Huber_Weight_calculation_IGRF14('SV')


%%% 0) General variables and paths
np = 10000; % number of spatial data points to evaluation at

radius = ones(np,1); %%% estimation is done at the normalized Earth's mean radius (Earth's surface), 
%%% change this for estimation at the CMBid points in space on which the weighting will be applied
rad = pi/180;

Pathout = '/users/yourname/work/IGRF14/Evaluation/Huber/'; 
%%% Pathout has to be changed for the outputs to be changed
%%% you should have this tree on you computer for the output figures
%%% '/users/yourname/work/IGRF14/Evaluation/Huber/SV/figures/'
%%% Similar directories have to be created for DGRF and IGRF

HOME = '/users/yourname/work/IGRF14/Evaluation/Huber/'; 
%%% Path HOME has to be change and be directed to the list of candidate models in the format
%%%[n m g h]
%%% SV models should be in HOME/SV/
%%% IGRF  should be in HOME/IGRF/
%%% DGRF should be in HOME/DGRF/

%%% CHECK inputs
if (TYPE == 'SV')
    PP = [HOME 'SV/'];
    Pathout = [Pathout '/SV/'];
    FILE = dir([PP '*.cof']);%%% scan available candidates
    Nmax = 8; % Maximum SH degree
    %%% LABELS have been manually defined by C. Beggan for IGRF-14
    LLf = {'BGS', 'CSES', 'DTU', 'Edinburgh', 'GFZ', 'IPGP', 'ISTERRE', 'Japan', ...
        'Leeds', 'MISTA', 'NASA', 'NOAA', 'OPGC', 'STRAS', 'TU_B', 'UCM', 'USTHB', 'WHU'};
    model_name = [Pathout 'SV_Huber.cof'];
    model_header = 'Candidate model with Huber weights in space, '; epoch = 2025;
elseif (TYPE == 'IG')
    PP = [HOME 'IGRF/'];
    Pathout = [Pathout 'IGRF/' ];
    FILE = dir([PP '*.cof']);%%% scan available candidates
    Nmax = 13; %%% Maximum SH degree for IGRF
    %%% LABELS have been defined by C. Beggan for IGRF-14
    LLf = {'BGS', 'CSES', 'DTU', 'GCRAS', 'GFZ', 'IPGP', 'ISTERRE', ... 
        'MISTA', 'NOAA', 'OPGC', 'STRAS', 'TU_B', 'UCM', 'USTHB', 'WHU'};
    model_name = [Pathout 'IGRF_Huber.cof'];
    model_header = 'Candidate model with Huber weights in space ';epoch = 2025;
elseif (TYPE == 'DG')
    PP = [HOME 'DGRF/'];
    Pathout = [Pathout 'DGRF/'];
    FILE = dir([PP '*.cof']);%%% scan available candidates
    Nmax = 13; %%% maximum SH degree for DGRF
    %%% LABELS have been defined by C. Beggan for IGRF-14
    LLf = {'BGS', 'CSES', 'DTU', 'GCRAS', 'GFZ', 'IPGP', 'ISTERRE', ... 
        'MISTA', 'NOAA', 'STRAS', 'TU_B', 'UCM', 'USTHB', 'WHU'};
    model_name = [Pathout 'DGRF_Huber.cof'];
    model_header = 'Candidate model with Huber weights in space ';epoch = 2020;
end

%%% 1) Generate an equal area grid with enough points to ensure
%%% orthogonality for the first maximum likelihood estimate
[Long Lat] = MakeAGrid(np);


TABCOEF = [];
for i=1:size(FILE,1),
    COEF = readmatrix([PP FILE(i).name], 'HeaderLines', 3, 'FileType', 'text');
    disp([PP FILE(i).name]);
    m = ConvertGHIntoVect(COEF);
    TABCOEF = [TABCOEF m];
    m_median = median(TABCOEF,2);
end

%%% 2) Compute the magnetic field of each candidate on the regular grid
%%%% 2-a) compute the design matrix
[A_r, A_th, A_phi] = design_SHA(radius, (90-Lat)*rad,Long*rad, Nmax);

Ball = [];
MatAll = [];

for i=1:size(TABCOEF,2)
    B_r(:,i) = A_r*TABCOEF(:,i);
    B_q(:,i) = A_th*TABCOEF(:,i);
    B_phi(:,i) = A_phi*TABCOEF(:,i);
    
    %%% Accumulate Matrix and data for Huber weighting
    Ball = [Ball;[B_r(:,i);B_q(:,i);B_phi(:,i)]];
    MatAll = [MatAll;[A_r;A_th;A_phi]];
end

%%% 2-b) Check that the conditionning of the design matrix is correct
%%% and does not introduce artificial numerical correlations between
%%% coefficients
MATA = MatAll'*MatAll;

%%% Preconditioning and normalization
DD=sqrt(diag(MATA)); V=diag(1./DD); U=diag(DD);
MKT=MATA*V; MK=V*MKT; clear MKT;
[~, S,~]=svd(MK);
LambMax=max(diag(S));
LambMin=min(diag(S));
clear S;
disp(['Regular Conditionning: ' num2str(cond(MATA))]);
disp(['Normalized Conditionning: ' num2str(LambMax/LambMin) ' (should be nearly equal to 1; if not then increase np).']);

%%% 3-c) Perform the IRLS with Huber weights in space
disp(['Performs Iteratively least-squares with Huber Weights in space using RobustFit.m']);

Nm = size(TABCOEF,2);
Grid = [Long Lat];
Nvd = size(Grid,1);

% USe the Matlab robustfit solveer
[m_i_huber,stats] =  robustfit(MatAll,Ball,'huber',1.345,'off');

%%% 3-d Get the weights after iteration for the plots
weights_mat = 1./(stats.w);

%%% 4) Prepare the grids of weights for each model
for i=1:Nm
    weights_huber{i} = [Grid weights_mat(1:Nvd) weights_mat(Nvd+1:2*Nvd) weights_mat(2*Nvd+1:3*Nvd)];
    weights_mat(1:3*Nvd) = [];
end

%%% 5) Plot the weigths and compute the RMS
% plot_weights(weights_huber,LLf,Pathout,'Weights'); % plot the figures if desired
TABLE = RMS_TABLE([TABCOEF median(TABCOEF,2) m_i_huber]); % compute RMS between the robust model, the candidates, and the median model
save([Pathout 'TABLE_RMS_with_Hubers.txt'],'TABLE','-ascii');
VOID = save_model_IGRF(m_i_huber, model_name,model_header,Nmax, epoch); % save the robust model
VOID = plot_res_map([m_i_huber m_median],13,'','MedianVsHuber',Pathout);
return



%%%%%%%%%%%%%%%%%%%%%%% ADDITIONAL EMBEDED FUNCTIONS
function VOID = Plot_mean_Weigth(GRID,weights_huber,Nm,Pathout);
VOID = [];
Weightsmean = zeros(length(GRID),1);

for i=1:Nm
    weightBr = weights_huber{i};
Weightsmean = (Weightsmean*(i-1)+1./weightBr(:,3))/i;
end
figure;
axesm mollweid;

% Define a colormap
scatterm(GRID(:,2),GRID(:,1),6,Weightsmean,'filled');
 map =[
        0.2510         0    0.2510
        0.2301         0    0.3134
        0.2092         0    0.3758
        0.1882         0    0.4383
        0.1673         0    0.5007
        0.1464         0    0.5631
        0.1255         0    0.6255
        0.1046         0    0.6879
        0.0837         0    0.7503
        0.0627         0    0.8127
        0.0418         0    0.8752
        0.0209         0    0.9376
        0         0    1.0000
        0.0111    0.1111    1.0000
        0.0222    0.2222    1.0000
        0.0333    0.3333    1.0000
        0.0444    0.4444    1.0000
        0.0556    0.5556    1.0000
        0.0667    0.6667    1.0000
        0.0778    0.7778    1.0000
        0.0889    0.8889    1.0000
        0.1000    1.0000    1.0000
        0.1900    1.0000    1.0000
        0.2800    1.0000    1.0000
        0.3700    1.0000    1.0000
        0.4600    1.0000    1.0000
        0.5500    1.0000    1.0000
        0.6400    1.0000    1.0000
        0.7300    1.0000    1.0000
        0.8200    1.0000    1.0000
        0.9100    1.0000    1.0000
        1.0000    1.0000    1.0000
        1.0000    1.0000    0.9000
        1.0000    1.0000    0.8000
        1.0000    1.0000    0.7000
        1.0000    1.0000    0.6000
        1.0000    1.0000    0.5000
        1.0000    1.0000    0.4000
        1.0000    1.0000    0.3000
        1.0000    1.0000    0.2000
        1.0000    1.0000    0.1000
        1.0000    1.0000         0
        1.0000    0.9091         0
        1.0000    0.8182         0
        1.0000    0.7273         0
        1.0000    0.6364         0
        1.0000    0.5455         0
        1.0000    0.4545         0
        1.0000    0.3636         0
        1.0000    0.2727         0
        1.0000    0.1818         0
        1.0000    0.0909         0
        1.0000         0         0
        1.0000         0    0.0909
        1.0000         0    0.1818
        1.0000         0    0.2727
        1.0000         0    0.3636
        1.0000         0    0.4545
        1.0000         0    0.5455
        1.0000         0    0.6364
        1.0000         0    0.7273
        1.0000         0    0.8182
        1.0000         0    0.9091
        1.0000         0    1.0000];
    
    colormap(map);
    colorbar
    caxis([0.5 1]);
    print('-djpeg','-r300', [Pathout 'Mean_HuberWeigths.jpg']);
    return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TABLE = RMS_TABLE(TAB)
for i=1:size(TAB,2)
    for j=1:size(TAB,2),
        S = powerspec(TAB(:,i)-TAB(:,j),1,'int'); TABLE(i,j) = sqrt(sum(S));
    end,
end
return



function VOID = plot_weights(TAB,LL,PATH,ARG)
%%% plot the weights for the three components on a global grid.
%%% define global grid
VOID = [];
C_max = 1;

load('~ciar/work/Swarm/coast.dat');
lat = coast(:,2);
long = coast(:,1);

for i=1:max(size(TAB))
    TAB_values = TAB{i};
   
   
clf;
for j=1:3
w_tmp = 1./TAB_values(:,j+2);
if (j==1)
    Comp = 'Br';
elseif(j==2)
    Comp = 'Btheta';
else
    Comp = 'Bphi';
end

 filename_out = [PATH ARG '_' LL{i} '_' Comp '.jpg'];
% North pole
pole_width = .3;
axes('Position', [0  0.45 pole_width pole_width], 'Box', 'off')
axesm('MapProjection', 'ortho', 'Origin', [90 0], 'FLatLimit',[-Inf 40.01])
plotm(lat, long, '-k', 'Linewidth', 1)
scatterm(TAB_values(:,2),TAB_values(:,1),15,w_tmp,'filled');
tightmap
caxis([-C_max C_max])
set(gca, 'Box', 'off', 'Visible', 'off')

axes('Position', [1-pole_width 0.45 pole_width pole_width ], 'Box', 'off')
axesm('MapProjection', 'ortho', 'Origin', [-90 0], 'FLatLimit',[-Inf 39.9])
plotm(lat, long, '-k', 'Linewidth', 1)
scatterm(TAB_values(:,2),TAB_values(:,1),15,w_tmp,'filled');
tightmap
caxis([-C_max C_max])
set(gca, 'Box', 'off', 'Visible', 'off')

axes('Position', [0 0 1 .5], 'Box', 'off')
axesm('MapProjection', 'hammer', 'Frame', 'off')
scatterm(TAB_values(:,2),TAB_values(:,1),15,w_tmp,'filled');
plotm(lat, long, '-k', 'Linewidth', 1)
gridm
tightmap
set(gca, 'Box', 'off', 'Visible', 'off')
caxis([-C_max C_max])
tightmap

% colorbar
axes('Position', [0.333 0.6 .333 .015], 'Box', 'off')

dC = 20;
set(gca,'XTick', [-C_max:dC:C_max], 'XtickLabel', [-C_max:dC:C_max], 'Xdir', 'Normal', 'FontWeight', 'normal', 'Xcolor', 'k', 'FontSize', 8, 'YTickLabel', {}, 'Ydir', 'Normal')
text(0.5, 2,'unit', 'Color', 'k', 'FontWeight','normal', 'FontSize', 10, 'Units', 'normalized', 'HorizontalAlignment', 'Center')


map =[
    0.2510         0    0.2510
    0.2301         0    0.3134
    0.2092         0    0.3758
    0.1882         0    0.4383
    0.1673         0    0.5007
    0.1464         0    0.5631
    0.1255         0    0.6255
    0.1046         0    0.6879
    0.0837         0    0.7503
    0.0627         0    0.8127
    0.0418         0    0.8752
    0.0209         0    0.9376
         0         0    1.0000
    0.0111    0.1111    1.0000
    0.0222    0.2222    1.0000
    0.0333    0.3333    1.0000
    0.0444    0.4444    1.0000
    0.0556    0.5556    1.0000
    0.0667    0.6667    1.0000
    0.0778    0.7778    1.0000
    0.0889    0.8889    1.0000
    0.1000    1.0000    1.0000
    0.1900    1.0000    1.0000
    0.2800    1.0000    1.0000
    0.3700    1.0000    1.0000
    0.4600    1.0000    1.0000
    0.5500    1.0000    1.0000
    0.6400    1.0000    1.0000
    0.7300    1.0000    1.0000
    0.8200    1.0000    1.0000
    0.9100    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    0.9000
    1.0000    1.0000    0.8000
    1.0000    1.0000    0.7000
    1.0000    1.0000    0.6000
    1.0000    1.0000    0.5000
    1.0000    1.0000    0.4000
    1.0000    1.0000    0.3000
    1.0000    1.0000    0.2000
    1.0000    1.0000    0.1000
    1.0000    1.0000         0
    1.0000    0.9091         0
    1.0000    0.8182         0
    1.0000    0.7273         0
    1.0000    0.6364         0
    1.0000    0.5455         0
    1.0000    0.4545         0
    1.0000    0.3636         0
    1.0000    0.2727         0
    1.0000    0.1818         0
    1.0000    0.0909         0
    1.0000         0         0
    1.0000         0    0.0909
    1.0000         0    0.1818
    1.0000         0    0.2727
    1.0000         0    0.3636
    1.0000         0    0.4545
    1.0000         0    0.5455
    1.0000         0    0.6364
    1.0000         0    0.7273
    1.0000         0    0.8182
    1.0000         0    0.9091
    1.0000         0    1.0000];
image([1:length(map)], 'XData', [0.2 C_max]);
colormap(map);
set(gcf, 'PaperOrientation', 'Portrait', 'PaperType', 'A4', 'PaperUnits', ...
    'centimeters', 'PaperPosition', [0.2 0.2 0.7*[30 30]]);

 print('-djpeg','-r300', filename_out);
end

end
return
 

function VOID = save_model_IGRF(m_i,model_name,model_header,maxdeg_MF,epoch)
% save model in IGRF format
N_koeff_MF  = maxdeg_MF*(maxdeg_MF+2);
fid=fopen(model_name,'w');
fprintf(fid,'# %s, %s\n', model_header, datestr(now));
fprintf(fid,'# Epoch = %7.2f  \n', epoch);
fprintf(fid,'# n   m        gnm        hnm       uncertainty_gnm     uncertainty_hnm\n');

i=1;
for n=1:maxdeg_MF;
    for m=0:n;
        if m == 0;
            fprintf(fid,'%3i %3i %10.2f %10.2f %10.2f %10.2f \n',n,m,m_i(i),0, 0.00, 0.00);
            i=i+1;
        elseif m > 0;
            
            fprintf(fid,'%3i %3i %10.2f %10.2f %10.2f %10.2f \n',n,m,m_i(i),m_i(i+1), 0.00, 0.00);
            
            i=i+2;
        end;
    end;
end;
fclose(fid);
VOID = [];
return

function m_i = read_model_igrf(filename, maxdeg_MF)
% [m_i, epoch, maxdeg_MF, maxdeg_sv, r_ref, m_ext] = read_model(filename)
% reads SHA coefficients stored in igrf format
%

if nargin < 2
if isempty(strfind(filename,'SV-2015-2020'))
  maxdeg_MF =13;
disp('Main field model')
else
    maxdeg_MF = 8;
    disp('Sv model');
end
end

fid=fopen(filename,'r');
header = fgetl(fid);
dummy = fgetl(fid);
dummy = fgetl(fid);
i=1;
N_koeff_MF = maxdeg_MF*(maxdeg_MF+2);
maxdeg_sv = 0;
for n=1:maxdeg_MF;
    for m=0:n;
        if n <= maxdeg_sv
            line_string = sscanf(fgetl(fid), '%i %i %f %f %f %f');
        else
            line_string = sscanf(fgetl(fid), '%i %i %f %f');
        end
        
        if (n ~= line_string(1)) & (m ~= line_string(2)) error 'n,m mismatch'; end;
        if m == 0;
            m_i(i) = line_string(3);
            if n <= maxdeg_sv
                m_i(N_koeff_MF+i) = line_string(5);
            end	
            i=i+1;
        elseif m > 0;
            m_i(i:i+1) = line_string(3:4);
            if n <= maxdeg_sv
                m_i(N_koeff_MF+i:N_koeff_MF+i+1) = line_string(5:6);
            end	
            i=i+2;
        end;
    end;
end;
m_ext = [];
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    m_ext = [m_ext; sscanf(tline, '%f', 1)];
end
fclose(fid);
return


function [Long,Lat]=MakeAGrid(N)
% [Long Lat]=MakeAGrid(N)
% Cree une grille reguliere sur la sphere en utilisant la routine matlab.
% Based on Ed Saff, "Sphere Points", 1997, uses partsphere.m
[Points,L,diam,topcol,botcol] = partsphere(N);

% Conversion de cartesien a spherique
P=Points';
[phi th r]=cart2sph(P(:,1),P(:,2),P(:,3));

% Conversion de radians a degres
Long=180*phi/pi;
Lat=180*th/pi;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% CODE FROM OTHER SOURCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Points,L,diam,topcol,botcol] = partsphere(N,lrounded,angles,mrounded,ipl)
%
% -- Input arguments --
%
% -- Output arguments --
%
%   References:
%
% [Points,L,diam] =  partsphere(N,lrounded,angles,mrounded,ipl)
%

%   $Revision: 1.5  $ Paul Leopardi 2003-10-13
%   Make angles=1 the default
%   $Revision: 1.4  $ Paul Leopardi 2003-09-18
%   Correct calculation of diameters
%   Return topcol and botcol: number of regions in top and bottom collars
%   $Revision: 1.3  $ Paul Leopardi 2003-09-13
%   Parameter angles controls whether to locate points via angle or height average
%   Match twist and offset to Maple version of the algorithm
%   $Revision: 1.2  $ Paul Leopardi 2003-09-09
%   Parameters lrounded and mrounded control rounding of L and m
%   $Revision: 1.1  $ Paul Leopardi 2003-09-08
%   Return diameters
%   $Revision: 1.0  $ Paul Leopardi 2003-09-08
%   for UNSW School of Mathematics

global area;
if nargin < 2
    lrounded = 0;
end
if nargin < 3
    angles = 1;
end
if nargin < 4
    mrounded = 0;
end
if nargin < 5
    ipl = 0;
end
if N == 1
    Points = zeros(3,N);
    Points(:,1)=[0,0,1]';
    L = 1;
    diam = 2*pi;
    topcol = 0;
    botcol = 0;
else
    area = 4*pi/N;
    % [A,val,exitflag,output] = fsolve(@square,sqrt(area),optimset(optimset('fsolve'),'Display','off'));
    % A = real(A);
    A = sqrt(area);
    Beta = acos(1-2/N);
    fuzz = eps*2*N;
    if lrounded == 1
        L = 2+max(1,round((pi-2*Beta)/A));
    else
        L = 2+max(1,ceil((pi-2*Beta)/A-fuzz));
    end
    Theta = (pi-2*Beta)/(L-2);
    mbar = zeros(1,L);
    mbar(1) = 1;
    for i = 2:L-1
        mbar(i) = N*(cos(Theta*(i-2)+Beta)-cos(Theta*(i-1)+Beta))/2;
    end
    mbar(L) = 1;
    alpha = 0;
    m = zeros(1,L);
    m(1) = 1;
    diam = zeros(1,L);
    diam(1) = 2*Beta;
    for i = 2:L
        if mrounded == 1
            m(i) = round(mbar(i)+alpha);
        else
            if (mbar(i)-floor(mbar(i)+alpha+fuzz)) < 0.5
                m(i) = floor(mbar(i)+alpha+fuzz);
            else
                m(i) = ceil(mbar(i)+alpha-fuzz);
            end
        end
        alpha = alpha + mbar(i)-m(i);
    end
    
    if ipl > 0
        fprintf('For N=%d, the number of levels cuts L=%d\n', N, L);
        fprintf('--------------------------------------------------\n');
    end
    Points = zeros(3,N);
    Points(:,1) = [0,0,1]';
    z = 1-(2+m(2))/N;
    Format = 'Level %3d: mbar(%2d) =%12.8f   m(%2d) =%4d diff=%12.8f; Area= %9.5f\n';
    i = 1;
    Area = 1;
    if ipl > 0
        fprintf(Format,i,i,mbar(i),i,m(i),mbar(i)-m(i),Area);
    end
    offset = zeros(1,L-1);
    offset(1) = 0;
    n=1;
    for i = 2:L-1
        twist=4;
        if m(i-1) ~= 0 & m(i) ~= 0
            offset(i) = offset(i-1)+gcd(m(i),m(i-1))/(2*m(i)*m(i-1))+min(twist,floor(m(i-1)/twist))/m(i-1);
        else
            offset(i) = 0;
        end
        top = z+m(i)/N;
        rtop = sqrt(1-top^2);
        bot = z-m(i)/N;
        rbot = sqrt(1-bot^2);
        if m(i) > 1
            angle = 2*pi/m(i);
            if rtop > rbot
                rside = rtop;
                hside = top;
            else
                rside = rbot;
                hside = bot;
            end
            side = acos([rside,0,hside]*[rside*cos(angle),rside*sin(angle),hside]');
            diag = acos([rtop,0,top]*[rbot*cos(angle),rbot*sin(angle),bot]');
            diam(i) = max(side,diag);
        else
            diam(i) = acos([rtop,0,top]*[-rbot,0,bot]');
        end
        if angles
            h = cos((acos(top)+acos(bot))/2);
        else
            h = z;
        end
        r = sqrt(1-h^2);
        for j = 0:m(i)-1
            s = offset(i)+j/m(i);
            n = n+1;
            Points(:,n) = [r*cos(2*pi*s),r*sin(2*pi*s),h]';
        end
        Area = 1;
        if ipl > 0
            fprintf(Format,i,i,mbar(i),i,m(i),mbar(i)-m(i),Area);
        end
        z = z-(m(i)+m(i+1))/N;
    end
    i = L;
    Area = 1;
    if ipl > 0
        fprintf(Format,i,i,mbar(i),i,m(i),mbar(i)-m(i),Area);
    end
    Points(:,N) = [0,0,-1]';
    diam(L) = 2*Beta;
    if L < 3
        topcol = 0;
        botcol = 0;
    else
        topcol = m(2);
        botcol = m(L-1);
    end
end
return

function v=square(A)
global area;
v=8*asin(sqrt(2)*sin(A/2)/sin(A))-2*pi-area;
return

function [A_r, A_theta, A_phi] = design_SHA(r, theta, phi, N, varargin);
% [A_r, A_theta, A_phi] = design_SHA(r, theta, phi, N)
%
% Calculates design matrices A_i that connects the vector
% of (Schmidt-normalized) spherical harmonic expansion coefficients,
% x = (g_1^0; g_1^1; h_1^1; g_2^0; g_2^1; h_2^1; ... g_N^N; h_N^N)
% and the magnetic component B_i, where "i" is "r", "theta" or "phi":
%        B_i = A_i*x
% Input: r(:)      radius vector (in units of the reference radius a)
%        theta(:)  colatitude    (in radians)
%        phi(:)    longitude     (in radians)
%        N         maximum degree/order
%
% [A_r, A_theta, A_phi] = design_SHA(r, theta, phi, N, i_e_flag)
% with i_e_flag = 'int' for internal sources (g_n^m and h_n^m)
%                 'ext' for external sources (q_n^m and s_n^m)
%
% Uses MEX file design_SHA_m if available, and Matlab program design_SHA_matlab.m else

% January 2003, Nils Olsen, DSRI
% Modified by E. Thebault for the time-varying version

% determine size of input arrays
max_size = max([size(r); size(theta); size(phi)]);
max_length = max_size(1)*max_size(2);
% convert to matrix if input parameter is scalar
if length(r)     == 1; r = r*ones(max_size); end;
if length(theta) == 1; theta = theta*ones(max_size); end;
if length(phi)   == 1; phi = phi*ones(max_size); end;
% check for equal length of all input parameters
if size(r) ~= size(theta) | size(r) ~= size(phi);
    error('Variables must be of equal size (or scalars)');
    return
end
% convert to row vector
r = reshape(r, max_length, 1);
theta = reshape(theta, max_length, 1);
phi = reshape(phi, max_length, 1);

if nargin == 4;
    i_e_flag = 'int'; % internal sources by default
elseif nargin > 4
    i_e_flag = varargin{1};
else
    error('At least 4 inputs required')
end;

[A_r, A_theta, A_phi] = design_SHA_matlab(r, theta, phi, N, i_e_flag);

return
% -----------------------------------------------------------------------

function [A_r, A_theta, A_phi] = design_SHA_matlab(r, theta, phi, N, varargin);
% [A_r, A_theta, A_phi] = design_SHA_matlab(r, theta, phi, N)
%
% Calculates design matrices A_i that connects the vector
% of (Schmidt-normalized) spherical harmonic expansion coefficients,
% x = (g_1^0; g_1^1; h_1^1; g_2^0; g_2^1; h_2^1; ... g_N^N; h_N^N)
% and the magnetic component B_i, where "i" is "r", "theta" or "phi":
%        B_i = A_i*x
% Input: r(:)      radius vector (in units of the reference radius a)
%        theta(:)  colatitude    (in radians)
%        phi(:)    longitude     (in radians)
%        N         maximum degree/order
%
% [A_r, A_theta, A_phi] = design_SHA(r, theta, phi, N, i_e_flag)
% with i_e_flag = 'int' for internal sources (g_n^m and h_n^m)
%                 'ext' for external sources (q_n^m and s_n^m)
%

% March 2001, Nils Olsen, DSRI

if nargin == 4;
    i_e_flag = 'int'; % internal sources by default
elseif nargin > 4
    i_e_flag = varargin{1};
else
    error('At least 4 inputs required')
end;
N_koeff=(N+1)^2-1;

cos_theta = cos(theta);
sin_theta = sin(theta);
LL = sin_theta == 0; sin_theta(LL) = 1e-5;
N_data    = length(theta);

A_r       = zeros(N_data, N_koeff);
A_theta   = zeros(N_data, N_koeff);
A_phi     = zeros(N_data, N_koeff);

k=0;
for n = 1:N
    if strcmp(i_e_flag, 'int')
        r_n       = r.^(-(n+2));
    elseif strcmp(i_e_flag, 'ext')
        r_n       = r.^(n-1);
    else
        warning 'i_e_flag neither "int" nor "ext". Assumed "int".'
        r_n       = r.^(-(n+2));
    end
    
    Pnm = legendre(n, cos_theta, 'sch')';      % P_n^m and derivatives vrt. theta
    dPnm(:,n+1) =  (sqrt(n/2).*Pnm(:,n));      % m=n
    dPnm(:,1) = -sqrt(n*(n+1)/2.).*Pnm(:,2);   % m=0
    if n > 1; dPnm(:,2)=(sqrt(2*(n+1)*n).*Pnm(:,1)-sqrt((n+2)*(n-1)).*Pnm(:,3))/2; end; % m=1
    for m = 2:n-1                              % m=2...n-1
        dPnm(:,m+1)=(sqrt((n+m)*(n-m+1)).*Pnm(:,m)-sqrt((n+m+1)*(n-m)).*Pnm(:,m+2))/2;
    end;
    if n == 1 dPnm(:,2) = sqrt(2)*dPnm(:,2); end
    
    if ~strcmp(i_e_flag, 'ext') % internal sources by default
        for m = 0:n;
            cos_phi   = cos(m*phi);
            sin_phi   = sin(m*phi);
            
            if m == 0
                k = k+1;       % terms corresponding to g_n^0
                A_r(:,k)       =  (n+1).*r_n(:).*Pnm(:,1);
                A_theta(:,k)   = -r_n(:).*dPnm(:,1);
                A_phi(:,k)     =  r_n(:)*0;
            else
                k = k+1;       % terms corresponding to g_n^m
                A_r(:,k)       =  (n+1).*r_n(:).*cos_phi.*Pnm(:,m+1);
                A_theta(:,k)   = -r_n(:).*cos_phi.*dPnm(:,m+1);
                A_phi(:,k)     =  r_n(:).*m.*sin_phi.*Pnm(:,m+1)./sin_theta;
                
                k = k+1;       % terms corresponding to h_n^m
                A_r(:,k)       =  (n+1).*r_n(:).*sin_phi.*Pnm(:,m+1);
                A_theta(:,k)   = -r_n(:).*sin_phi.*dPnm(:,m+1);
                A_phi(:,k)     = -r_n(:).*m.*cos_phi.*Pnm(:,m+1)./sin_theta;
            end
        end; % m
    else % external sources, i_e_flag ='ext'
        for m = 0:n;
            cos_phi   = cos(m*phi);
            sin_phi   = sin(m*phi);
            
            if m == 0
                k = k+1;       % terms corresponding to q_n^0
                A_r(:,k)       = -n.*r_n(:).*Pnm(:,1);
                A_theta(:,k)   = -r_n(:).*dPnm(:,1);
                A_phi(:,k)     =  r_n(:)*0;
            else
                k = k+1;       % terms corresponding to q_n^m
                A_r(:,k)       = -n.*r_n(:).*cos_phi.*Pnm(:,m+1);
                A_theta(:,k)   = -r_n(:).*cos_phi.*dPnm(:,m+1);
                A_phi(:,k)     =  r_n(:).*m.*sin_phi.*Pnm(:,m+1)./sin_theta;
                
                k = k+1;       % terms corresponding to s_n^m
                A_r(:,k)       = -n.*r_n(:).*sin_phi.*Pnm(:,m+1);
                A_theta(:,k)   = -r_n(:).*sin_phi.*dPnm(:,m+1);
                A_phi(:,k)     = -r_n(:).*m.*cos_phi.*Pnm(:,m+1)./sin_theta;
            end
        end; % m
    end; % if int/ext
end; % n

return

function VOID = plot_res_map(TAB,N,LL,ARG,PATH,LEG)

%%% define global grid
Step = 2;
VOID = [];
theta = [180-Step/2:-Step:0];
phi   = [-180+Step/2:Step:180-Step/2];

a = 6371.2;
C_max = 10; dC = 2;
load('~ciar/work/Swarm/coast.dat')

lat = coast(:,2);
long = coast(:,1);

if size(TAB,2) ~= 2
    for i=1:size(TAB,2)
        mi = TAB(:,i);
        
        dB_r = synth_grid(mi,(a+00)/a, theta, phi);
        %dB_r = atan2(Bp,-Bq)*180/pi;
        
        %B_ref =  makerefmat('RasterSize', size(dB_r), ...
        %    'Latlim', [-90 90], 'Lonlim', [-180 180]);
        B_ref = georefcells([-90 90],[-180 180], size(dB_r));
        
        clf;
        B_tmp = dB_r;
        filename_out = [PATH 'dBr_' LL{i} '_with_' ARG '.jpg'];
        % North pole
        pole_width = .3;
        axes('Position', [0  0.45 pole_width pole_width], 'Box', 'off')
        axesm('MapProjection', 'ortho', 'Origin', [90 0], 'FLatLimit',[-Inf 40.01])
        plotm(lat, long, '-k', 'Linewidth', 1)
        meshm(B_tmp, B_ref)
        tightmap
        caxis(C_max*[-1 1])
        set(gca, 'Box', 'off', 'Visible', 'off')
        
        axes('Position', [1-pole_width 0.45 pole_width pole_width ], 'Box', 'off')
        axesm('MapProjection', 'ortho', 'Origin', [-90 0], 'FLatLimit',[-Inf 39.9])
        plotm(lat, long, '-k', 'Linewidth', 1)
        meshm(B_tmp, B_ref)
        tightmap
        caxis(C_max*[-1 1])
        set(gca, 'Box', 'off', 'Visible', 'off')
        
        axes('Position', [0 0 1 .5], 'Box', 'off')
        axesm('MapProjection', 'hammer', 'Frame', 'off')
        meshm(B_tmp, B_ref)
        plotm(lat, long, '-k', 'Linewidth', 1)
        gridm
        tightmap
        set(gca, 'Box', 'off', 'Visible', 'off')
        caxis(C_max*[-1 1])
        tightmap
        
        % colorbar
        axes('Position', [0.333 0.6 .333 .015], 'Box', 'off')
        
        
        map =[
            0.2510         0    0.2510
            0.2301         0    0.3134
            0.2092         0    0.3758
            0.1882         0    0.4383
            0.1673         0    0.5007
            0.1464         0    0.5631
            0.1255         0    0.6255
            0.1046         0    0.6879
            0.0837         0    0.7503
            0.0627         0    0.8127
            0.0418         0    0.8752
            0.0209         0    0.9376
            0         0    1.0000
            0.0111    0.1111    1.0000
            0.0222    0.2222    1.0000
            0.0333    0.3333    1.0000
            0.0444    0.4444    1.0000
            0.0556    0.5556    1.0000
            0.0667    0.6667    1.0000
            0.0778    0.7778    1.0000
            0.0889    0.8889    1.0000
            0.1000    1.0000    1.0000
            0.1900    1.0000    1.0000
            0.2800    1.0000    1.0000
            0.3700    1.0000    1.0000
            0.4600    1.0000    1.0000
            0.5500    1.0000    1.0000
            0.6400    1.0000    1.0000
            0.7300    1.0000    1.0000
            0.8200    1.0000    1.0000
            0.9100    1.0000    1.0000
            1.0000    1.0000    1.0000
            1.0000    1.0000    0.9000
            1.0000    1.0000    0.8000
            1.0000    1.0000    0.7000
            1.0000    1.0000    0.6000
            1.0000    1.0000    0.5000
            1.0000    1.0000    0.4000
            1.0000    1.0000    0.3000
            1.0000    1.0000    0.2000
            1.0000    1.0000    0.1000
            1.0000    1.0000         0
            1.0000    0.9091         0
            1.0000    0.8182         0
            1.0000    0.7273         0
            1.0000    0.6364         0
            1.0000    0.5455         0
            1.0000    0.4545         0
            1.0000    0.3636         0
            1.0000    0.2727         0
            1.0000    0.1818         0
            1.0000    0.0909         0
            1.0000         0         0
            1.0000         0    0.0909
            1.0000         0    0.1818
            1.0000         0    0.2727
            1.0000         0    0.3636
            1.0000         0    0.4545
            1.0000         0    0.5455
            1.0000         0    0.6364
            1.0000         0    0.7273
            1.0000         0    0.8182
            1.0000         0    0.9091
            1.0000         0    1.0000];
        
        colormap(map);
        
        image([1:length(map)], 'XData', [-C_max C_max]);
        set(gca,'XTick', [-C_max:dC:C_max], 'XtickLabel', [-C_max:dC:C_max], 'Xdir', 'Normal', 'FontWeight', 'normal', 'Xcolor', 'k', 'FontSize', 8, 'YTickLabel', {}, 'Ydir', 'Normal')
        %text(0.5, 2,LEG, 'Color', 'k', 'FontWeight','normal', 'FontSize', 10, 'Units', 'normalized', 'HorizontalAlignment', 'Center')
        set(gcf, 'PaperOrientation', 'Portrait', 'PaperType', 'A4', 'PaperUnits', ...
            'centimeters', 'PaperPosition', [0.2 0.2 0.7*[30 30]]);
        
        print('-djpeg','-r300', filename_out);
    end
    
else
    depth = 0;
    [dB_r1,Bq1,Bp1] = synth_grid(TAB(:,1),(a+depth)/a, theta, phi);
    D1 = atan2(Bp1,-Bq1)*180/pi;
    
    [dB_r2,Bq2,Bp2] = synth_grid(TAB(:,2),(a+depth)/a, theta, phi);
    D2 = atan2(Bp2,-Bq2)*180/pi;
    %dB_r = (D1-D2)*180/pi;
    dB_r = dB_r2-dB_r1;
    %B_ref =  makerefmat('RasterSize', size(dB_r), ...
     %   'Latlim', [-90 90], 'Lonlim', [-180 180]);
     B_ref = georefcells([-90 90],[-180 180], size(dB_r));


    clf;
    B_tmp = dB_r;
    filename_out = [PATH 'Br_with_' ARG '.jpg'];
    % North pole
    pole_width = .3;
    axes('Position', [0  0.45 pole_width pole_width], 'Box', 'off')
    axesm('MapProjection', 'ortho', 'Origin', [90 0], 'FLatLimit',[-Inf 40.01])
    plotm(lat, long, '-k', 'Linewidth', 1)
    meshm(B_tmp, B_ref)
    tightmap
    caxis(C_max*[-1 1])
    set(gca, 'Box', 'off', 'Visible', 'off')
    
    axes('Position', [1-pole_width 0.45 pole_width pole_width ], 'Box', 'off')
    axesm('MapProjection', 'ortho', 'Origin', [-90 0], 'FLatLimit',[-Inf 39.9])
    plotm(lat, long, '-k', 'Linewidth', 1)
    meshm(B_tmp, B_ref)
    tightmap
    caxis(C_max*[-1 1])
    set(gca, 'Box', 'off', 'Visible', 'off')
    
    axes('Position', [0 0 1 .5], 'Box', 'off')
    axesm('MapProjection', 'hammer', 'Frame', 'off')
    meshm(B_tmp, B_ref)
    plotm(lat, long, '-k', 'Linewidth', 1)
    gridm
    tightmap
    set(gca, 'Box', 'off', 'Visible', 'off')
    caxis(C_max*[-1 1])
    tightmap
    
    % colorbar
    axes('Position', [0.333 0.6 .333 .015], 'Box', 'off')
    
    
    map =[
        0.2510         0    0.2510
        0.2301         0    0.3134
        0.2092         0    0.3758
        0.1882         0    0.4383
        0.1673         0    0.5007
        0.1464         0    0.5631
        0.1255         0    0.6255
        0.1046         0    0.6879
        0.0837         0    0.7503
        0.0627         0    0.8127
        0.0418         0    0.8752
        0.0209         0    0.9376
        0         0    1.0000
        0.0111    0.1111    1.0000
        0.0222    0.2222    1.0000
        0.0333    0.3333    1.0000
        0.0444    0.4444    1.0000
        0.0556    0.5556    1.0000
        0.0667    0.6667    1.0000
        0.0778    0.7778    1.0000
        0.0889    0.8889    1.0000
        0.1000    1.0000    1.0000
        0.1900    1.0000    1.0000
        0.2800    1.0000    1.0000
        0.3700    1.0000    1.0000
        0.4600    1.0000    1.0000
        0.5500    1.0000    1.0000
        0.6400    1.0000    1.0000
        0.7300    1.0000    1.0000
        0.8200    1.0000    1.0000
        0.9100    1.0000    1.0000
        1.0000    1.0000    1.0000
        1.0000    1.0000    0.9000
        1.0000    1.0000    0.8000
        1.0000    1.0000    0.7000
        1.0000    1.0000    0.6000
        1.0000    1.0000    0.5000
        1.0000    1.0000    0.4000
        1.0000    1.0000    0.3000
        1.0000    1.0000    0.2000
        1.0000    1.0000    0.1000
        1.0000    1.0000         0
        1.0000    0.9091         0
        1.0000    0.8182         0
        1.0000    0.7273         0
        1.0000    0.6364         0
        1.0000    0.5455         0
        1.0000    0.4545         0
        1.0000    0.3636         0
        1.0000    0.2727         0
        1.0000    0.1818         0
        1.0000    0.0909         0
        1.0000         0         0
        1.0000         0    0.0909
        1.0000         0    0.1818
        1.0000         0    0.2727
        1.0000         0    0.3636
        1.0000         0    0.4545
        1.0000         0    0.5455
        1.0000         0    0.6364
        1.0000         0    0.7273
        1.0000         0    0.8182
        1.0000         0    0.9091
        1.0000         0    1.0000];
    
    colormap(map);
  
    image([1:length(map)], 'XData', [-C_max C_max]);
    set(gca,'XTick', [-C_max:dC:C_max], 'XtickLabel', [-C_max:dC:C_max], 'Xdir', 'Normal', 'FontWeight', 'normal', 'Xcolor', 'k', 'FontSize', 8, 'YTickLabel', {}, 'Ydir', 'Normal')
    %text(0.5, 2,LEG, 'Color', 'k', 'FontWeight','normal', 'FontSize', 10, 'Units', 'normalized', 'HorizontalAlignment', 'Center')
    set(gcf, 'PaperOrientation', 'Portrait', 'PaperType', 'A4', 'PaperUnits', ...
        'centimeters', 'PaperPosition', [0.2 0.2 0.7*[30 30]]);
      title(ARG);
       print('-djpeg','-r300', filename_out);
    
end
return
 


function gh = ConvertGHIntoVect(GH)
gh = [];
for n=1: max(GH(:,1))
    for m=0:n
    if (m==0)
        gh = [gh;GH(1,3)];
        GH(1,:)=[];
    else
    gh = [gh;GH(1,3:4)'];
    GH(1,:) = [];
    end
    
    end
end
return

function R_n = powerspec(m_i, varargin)
% R_n = powerspec(m_i)
% Lowe-Mauersberger spectrum of Gauss coefficients m_i at radius r/a = 1
% R_n = powerspec(m_i, rho, 'int') yields spectrum at radius r/a = rho
% 
% May 2000, Nils Olsen, DSRI, Copenhagen

if nargin == 1
   rho = 1;
else
   rho = varargin{1};
end

if nargin > 1
    source = varargin{2};
else
    source = 'int'; % internal sources by default
end
    
maxdeg_MF = sqrt(length(m_i)+1)-1;
R_n = zeros(maxdeg_MF, 1);
i=0;
for n = 1:maxdeg_MF
   for m = 0:n
      i=i+1;
      R_n(n) = R_n(n)+m_i(i).^2;
      if m > 0
         i=i+1;
         R_n(n) = R_n(n)+m_i(i).^2;
      end
   end
   if strcmp(source, 'int')   
       R_n(n) = R_n(n)*(n+1)*rho.^-(2*n+4);
   elseif strcmp(source, 'ext')   
       R_n(n) = R_n(n)*n*rho.^(2*n-4);
   elseif strcmp(source, 'tor')   
       R_n(n) = R_n(n)*n*(n+1)/(2*n+1)*rho.^-2;
   else
       warning 'Unrecognized source: not ''int'', ''ext'', or ''tor'''
   end
end
function [B_r, varargout] = synth_grid(m_i, r, theta, phi, varargin)
% B_r                      = synth_grid(m_i, r, theta, phi)
% [B_r, B_theta]           = synth_grid(m_i, r, theta, phi)
% [B_r, B_theta, B_phi]    = synth_grid(m_i, r, theta, phi)
% [B_r, B_theta, B_phi, F] = synth_grid(m_i, r, theta, phi)
% [B_r, B_theta, B_phi, F] = synth_grid(m_i, r, theta, phi, source)
% 
% Synthetic magnetic field values B_r, B_theta, B_phi, F (in nT) 
% on a regular grid given by radius r (in units of the Earth's radius),
% co-latitude theta(:) and longitude phi(:) (in degrees)
% m_i(:) contains the Gauss coefficients (in nT, Schmidt normalized)
% source: optional 5th argrument specifies if coefficients are for
% internal ('int', default), external ('ext') or toroidal ('tor') sources
%
% July 2003, Nils Olsen, DSRI

% Last change: NIO 080703: only the specified field components are calculated

theta = theta(:)';
phi = phi(:)';

rad = pi/180;
maxdeg_MF = sqrt(length(m_i)+1)-1;
cos_theta = cos(theta*rad);
sin_theta = sin(theta*rad);
N_theta   = length(theta);
N_phi     = length(phi);
n = [1:maxdeg_MF];

if nargin > 4
    source = varargin{1};
else
    source = 'int'; % internal sources by default
end
if strcmp(source, 'int')   
    r_n       = r.^(-(n+2));
    fac_r     = n+1;
elseif strcmp(source, 'ext')   
    r_n       = r.^(+(n-1));
    fac_r     = -n;
elseif strcmp(source, 'tor')   
    r_n       = repmat(r.^(1), size(n));
    fac_r     = -n.*(n+1);
else
    warning 'Unrecognized source: not ''int'', ''ext'', or ''tor'''
end

n_koeff = maxdeg_MF*(maxdeg_MF+2);
T_r     = zeros(N_theta, n_koeff);
cos_sin_m = zeros(N_phi, n_koeff);
if nargout > 1; % also calculation of B_theta
    T_theta = zeros(N_theta, n_koeff);
end;
if nargout > 2; % also calculation of B_phi
    T_phi   = zeros(N_theta, n_koeff);
    sin_cos_m = zeros(N_phi, n_koeff);
end

k=1;
for n = 1:maxdeg_MF  
    Pnm = legendre(n, cos_theta, 'sch')';
    
    if nargout > 1; % also calculation of B_theta
        dPnm(:,n+1) =  (sqrt(n/2).*Pnm(:,n));      % m=n
        dPnm(:,1) = -sqrt(n*(n+1)/2.).*Pnm(:,2);   % m=0
        if n > 1; dPnm(:,2)=(sqrt(2*(n+1)*n).*Pnm(:,1)-sqrt((n+2)*(n-1)).*Pnm(:,3))/2; end; % m=1
        for m = 2:n-1                              % m=2...n-1
            dPnm(:,m+1)=(sqrt((n+m)*(n-m+1)).*Pnm(:,m)-sqrt((n+m+1)*(n-m)).*Pnm(:,m+2))/2;
        end;
        if n == 1 dPnm(:,2) = sqrt(2)*dPnm(:,2); end
    end
    
    T_r(:,k)     = fac_r(n).*r_n(n).*Pnm(:,1);
    cos_sin_m(:,k) = ones(N_phi, 1);

    if nargout > 1; % also calculation of B_theta
        T_theta(:,k) = -r_n(n).*dPnm(:,1);
    end
    if nargout > 2; % also calculation of B_phi
        T_phi(:,k)   = zeros(N_theta, 1);
        sin_cos_m(:,k) = zeros(N_phi, 1);
    end
    k=k+1;
    for m = 1:n
        T_r(:,k:k+1)     =  fac_r(n).*r_n(n).*[Pnm(:,m+1) Pnm(:,m+1)];
        cos_sin_m(:,k:k+1) = [cos(m*phi'*rad)  sin(m*phi'*rad)];
        if nargout > 1; % also calculation of B_theta
            T_theta(:,k:k+1) = -r_n(n).*[dPnm(:,m+1) dPnm(:,m+1)];
        end
        if nargout > 2; % also calculation of B_phi
            T_phi(:,k:k+1)   = m.*r_n(n).*[Pnm(:,m+1)./sin_theta' Pnm(:,m+1)./sin_theta'];
            sin_cos_m(:,k:k+1) = [sin(m*phi'*rad) -cos(m*phi'*rad)];
        end
        k = k+2;
    end
end

m = spdiags(m_i, 0, n_koeff, n_koeff);
B_r = T_r*(cos_sin_m*m)';

if nargout > 1; % also calculation of B_theta
    varargout{1} = T_theta*(cos_sin_m*m)';
end;

if nargout > 2; % also calculation of B_phi
    varargout{2} = T_phi*(sin_cos_m*m)';
end;

if nargout == 4;  % also calculation of F
    varargout{3} = sqrt(B_r.^2 + varargout{1}.^2 + varargout{2}.^2); 
end;
