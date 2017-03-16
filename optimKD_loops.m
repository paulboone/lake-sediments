
%% ttlem optimization procedure of K and D based on a single lake
%% procedure to run model for optimization:
%>compute the observed erosion rate  and its uncertainty
%>define the extent of the basin area as well as the lake elements
% >make an array of K D combinations.
% >run the model to a recent amount of time where I have volume of sed
% computed
% > compute the volume of seds produced - the total change in elevation
% within basin and not in lake*area
% >compute the chi^2 between modeled and observed, for the error, use the
% v/t err. Which is: 
% the 

%%
%% >load the DEM
%{
PREREQ to even loading the file:
Note that, throughout the use of TopoToolbox, it is assumed that the DEM has 
a projected coordinate system (e.g. UTM WGS84) and that elevation and horizontal 
coordinates are in meter units.
%}
close all, clear all, clc % clear all figures, all variables, and command window
fpath='/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/';
fname='clipped_jan1.tif';
name=[fpath,fname]; % char concatenation
DEM = GRIDobj(name);



%% >compute the observed erosion rate  and its uncertainty
t0=6000%yr
t1=0
D=1.73;%m deposits

%> DEM smooth
hsize=10;%15
sigma=2;%2
h = fspecial('gaussian', hsize, sigma); % 2D gaussian filter from image processing toolkit
%surf(h),shg

DEM2=DEM;
DEM2.Z=imfilter(DEM.Z, h, 'replicate'); % applies the image filter
DEM=DEM2;

%%fill sinks to correct for 'erroneous topographic depressions'
DEMf = fillsinks(DEM); % topotoolbox

%%flow direction object:
FD = FLOWobj(DEMf); % topotoolbox

%%flow acc
A = flowacc(FD);  % topotoolbox
G=gradient8(DEM); % topotoolbox
    
imageschs(DEM, min(gradient8(DEM),1)) % image scaled colors  & hillshade

%%when outlet is already known

% imageschs(DEM, min(gradient8(DEM),1)
i=15;%jan
load ('resM_pairs19')
resM=resM_pairs19;

x=resM(i,1)
y=resM(i,2)

%% compute drainage basin for this outlet
IX=coord2ind(DEM,x,y);      % convert xy coordinates (NOT indices) to linear index to get Z index of outlet
zo=DEM.Z(IX);               % outlet elev
B=drainagebasins(FD,IX);    % get a binary grid of basin cells;
   
%% detect lake 
%{
create binary grid of lake with the conditions:
- must be in drainage basin
- must have gradient < 0.03 <= WHAT does this actually mean? is this saying
it is ~ flat?
- must have height < less than outlet height + 1 (meter??)
%}
lk=DEM.Z*0;
lk(B.Z==1 & G.Z<0.03 & DEM.Z<zo+1)=1; % do all these conditions have a specific source?

% create grid of z heights for every point in drainage basin
Zb=DEM.Z*NaN;
Zb(B.Z==1)=DEM.Z(B.Z==1);

% plot z heights for drainage basin
imagesc(Zb);
axis equal
shg

%% vol of sed in time frame
Al=length(find(lk==1)); % # of grid cells that are in the lake
Aw=A.Z(IX)-Al; % # of upstream grid cells from outlet - cells in lake
% shouldn't this be the same size as the basin though??

pa=DEM.cellsize*DEM.cellsize; % cell area
Vemp=D*Al*pa;   % total sediment volume = meters of sediment deposits * cells in lake * area in cell
Vemp_u=Vemp/2;  % the case of triangular cs rather than rect ??

%% >make an array of K D combinations.

%%Temporal domain for real run
p.TimeSpan=t0;
p.TimeStep=t0;

% %%temporal domain for example run
% p.TimeSpan=1e6;
% p.TimeStep=1e5;

%%Hillslope processes, parameters values
%Set diffusivity parameter
p.D=0.03;

%%River incision parameters
p.Kw = 3e-6;
p.m  = 0.5;
p.n  = 1;
p.AreaThresh = 2e5; % channel contributing area threshold, m²

% (un)comment to get the drainage development type of choice
% Fixed or variable drainage netwrok through time
 p.DrainDir='variable';
% p.DrainDir='fixed';


%%Threshold slopes
p.Sc=1; %%21°) Not the default, why?
p.Sc_unit='tangent';


%%Output
p.ploteach=1000000;
p.saveeach=1;

% Specify the location where the results can be stored (e.g.
% p.resultsdir='C:\...\'); The default folder is the result file where the
% main model structure is stored.
% p.resultsdir='C:\...';
p.fileprefix=['JanRunK',num2str(p.Kw), 'D', num2str(p.D), 't'];

%% Initialize parameter structure.
% By making p an instance of ttlemset, the user ensures parameter values
% are set in the right way
p   = ttlemset(p);
%% >run the model to a recent amount of time where I have volume of sed
% computed
U=DEM;
U.Z=DEM.Z*0;
ttlem_out = ttlem(DEM,U,p);

%% > compute the volume of seds produced - the total change in elevation
% within basin and not in lake*area
out_name=[p.resultsdir 'JanRunK',num2str(p.Kw), 'D', num2str(p.D),'t', num2str(t0), '.mat'];
load(out_name,'H1');
h0=H1;

%%
close all
dzdt=DEM.Z-h0.Z; % difference between original and evolved DEM

dzdtb=dzdt*NaN;
dzdtb(B.Z==1 & lk~=1) = dzdt(B.Z==1 & lk~=1); % create ?z grid for cells in basin but not in lake

subplot(1,2,1)
imagesc(DEM),shg % show original DEM
subplot(1,2,2)   
imagesc(dzdtb),shg % show ?z grid DEM
axis image;

%% > compute the total eroded volume
Vmod=nansum(dzdtb(:))*pa


% >compute the chi^2 between modeled and observed, for the error, use the
% v/t err. Which is: 
% the 

chi2=(Vmod-Vemp).^2/Vemp_u^2;


%% now run optim loop to ehceck procedure

% D=[0:0.05:0.25]*10
% K=[3e-5:0.0005:3e-3]*10
D=linspace(0, 0.5,20);      %0:0.05:2.5*10
K=linspace(0, 0.025,20);    %3e-5:0.0005:3e-2


%
Mk=0;
Md=0;
chi2=0;
for i=1:length(D);
    for j=1:length(K)
%         111
        Mk(i,j)=K(j);
        Md(i,j)=D(i);
        %hillslope diff
        p.D=D(i);
        %%River incision parameters
        p.Kw =K(j);
%         222
        %%
        p.fileprefix=['JanRunK',num2str(p.Kw), 'D', num2str(p.D), 't']; ;
        
        %% Initialize parameter structure.
        % By making p an instance of ttlemset, the user ensures parameter values
        % are set in the right way
        p   = ttlemset(p);
        %% >run the model to a recent amount of time where i have volume of sed
        % computed
        U=DEM;
        U.Z=DEM.Z*0;
        ttlem_out = ttlem(DEM,U,p);
%         333
        %% > compute the volume of seds produced - the total change in elevation
        % within basin and not in lake*area
        out_name=[p.resultsdir 'JanRunK',num2str(p.Kw), 'D', num2str(p.D),'t', num2str(t0), '.mat'];
        load(out_name,'H1');
        h0=H1;
        
        %%
        close all
        dzdt=DEM.Z-h0.Z;
        dzdtb=dzdt*NaN;
        dzdtb(B.Z==1 & lk~=1)=dzdt(B.Z==1 & lk~=1);
        
%         subplot(1,2,1)
%         imagesc(DEM),shg
%         subplot(1,2,2)
%         imagesc(dzdtb),shg
%         axis image;
        
        %% > compute the eroded volume
        Vmod=nansum(dzdtb(:))*pa
        
        
        % >compute the chi^2 between modeled and observed, for the error, use the
        % v/t err. Which is:
        % the
        
        chi2(i,j)=(Vmod-Vemp).^2/Vemp_u^2;
        
        disp(j)
    end
    disp('#####################################')
    disp(i)
    disp('#####################################')
end

%%
% close all, contour(log10(chi2),20),shg  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ttlem optimization procedure of K and D based on 2 lakes lake
%% procedure to run model for optimization:
%>compute the observed erosion rate  and its uncertainty
%>define the extent of the basin area as well as the lake elements
% >make an array of K D combinations.
% >run the model to a recent amount of time where I have volume of sed
% computed
% > compute the volume of seds produced - the total change in elevation
% within basin and not in lake*area
% >compute the chi^2 between modeled and observed, for the error, use the
% v/t err. Which is: 
% the 

%% make dem for lost lake
%%>load the DEM 1
close all, clear all, clc
fpath='/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/';
fname='clipped_jan1.tif';
name=[fpath,fname];
DEM1 = GRIDobj(name);

%%get outlet 1
%%when outlet is already known
i=15;%jan
load ('resM_pairs19')
resM=resM_pairs19;

x1=resM(i,1);
y1=resM(i,2);


%% >load the DEM 2
fpath='/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/';
fname='clipped_lost.tif';
name=[fpath,fname];
DEM2 = GRIDobj(name);


% close all
% 
%     DEMs=DEM2;
%     DEMs.Z=imfilter(DEM.Z, h, 'replicate');
%     DEM=DEMs;
%     
%     %%fill sinks
%     DEMf = fillsinks(DEM);
%     
%     %%flow direction object:
%     FD = FLOWobj(DEMf);
%     
%     %%flow acc
%     A = flowacc(FD);
%     imageschs(DEM,dilate(A.^0.3,ones(1)),'colormap',flipud(copper));
%     hold on
%     contour(DEM, 60, '-k')
%     shg;
%     pause%this is for zooming in
%     %% get outlet after zooming
%     [x3,y3]=ginput(1)

%%
x2=3.523578530974025e+05;
y2=1.612583134776515e+06;

% Dv=[D1, D2];

%% inter the thickness at time at lost approxi 138
% close all
% plot ([122, 158], 	[5395,6607])
% shg

%% >compute the observed erosion rate  and its uncertainty

close all
t0=6000%yr
t1=0
D1=1.73;%jan m deposits
D2=1.38;%lost m deposits - interp
Dv=[D1, D2];

DEMv={DEM1, DEM2};
xv=[x1, x2];
yv=[y1, y2];
Bv={0};
lkv={0};
%% compute erosion volume
% prepare DEMs
for i=1:2%loop through to compute erosion volume from real data
    %loop between dems to opimize chi 2 by 2 dems%the optim K for jan is 0.001, teh optim D is 0.2
    %> DEM smooth
    hsize=10;%15
    sigma=2;%2
    h = fspecial('gaussian', hsize, sigma);
    %surf(h),shg
    1111
    
    DEM=DEMv{i};
    DEMs=DEM;
    DEMs.Z=imfilter(DEM.Z, h, 'replicate');
    DEM=DEMs;
    
    %%fill sinks
    DEMf = fillsinks(DEM);
    
    %%flow direction object:
    FD = FLOWobj(DEMf);
    
    %%flow acc
    A = flowacc(FD);
    G=gradient8(DEM);
    
    2222
    %%compute drainage basin for this outlet
    x=xv(i);
    y=yv(i);
    IX=coord2ind(DEM,x,y);%get index of outlet
    
    zo=DEM.Z(IX);%outlet elev
    B=drainagebasins(FD,IX);%get a binary grid of basin cells;
    Bv{i}=B;
    %%detect lake
    lk=DEM.Z*0;
    lk(B.Z==1 & G.Z<0.03 & DEM.Z<zo+1)=1;
    lkv{i}=lk;
    %
    Zb=DEM.Z*NaN;
    Zb(B.Z==1)=DEM.Z(B.Z==1);
    %
    imagesc(Zb);
    axis equal
    shg
%     pause
    
    3333
    %%vol of sed in time frame
    Al=length(find(lk==1));
    Aw=A.Z(IX)-Al
    Aw=nansum(B.Z(B.Z==1 & lk~=1))
    pa=DEM.cellsize*DEM.cellsize;
    %%
    D=Dv(i);
    Vemp(i)=D*Al*pa;
    Vemp_u(i)=Vemp(i)/2; %teh case of traingular cs rather than rect
end


%% >make an array of K D combinations.

%%Temporal domain
p.TimeSpan=t0;
p.TimeStep=t0;



%%Hillslope processes, parameters values
p.m  = 0.5;
p.n  = 1;

% (un)comment to get the drainage development type of choice
% Fixed or variable drainage netwrok through time
 p.DrainDir='variable';
% p.DrainDir='fixed';


%%Threshold slopes
p.Sc=1; %%21°)
p.Sc_unit='tangent';


%%Output
p.ploteach=inf;
p.saveeach=1;



%% now run optim loop to check procedure

% D=[0:0.05:0.25]*10
% K=[3e-5:0.0005:3e-3]*10
%% with A thresh
p.AreaThresh = 2e5; % channel contributing area threshold, m²
D=linspace(0,0.3,10);%0:0.05:2.5*10
K=linspace(0, 0.01,10);%3e-5:0.0005:3e-2

% %% without A thresh
% close all
% p.AreaThresh = 0; % channel contributing area threshold, m²
% % D=linspace(0.3,1.2,10);%0:0.05:2.5*10
% % K=linspace(0, 0.00001,10);%3e-5:0.0005:3e-2
% D=linspace(0.9,1.2,10);%0:0.05:2.5*10
% K=linspace(0, 0.0000002,10);%3e-5:0.0005:3e-2



[Md, Mk]=meshgrid(D,K);

%
Mk=0;
Md=0;
chi2M=0;

namev={'jan', 'lost'};

%jan is 0.001, teh optim D is 0.2

for i=1:length(D);
    for j=1:length(K)
        
        111
        Mk(i,j)=K(j);
        Md(i,j)=D(i);
        %hillslope diff
        p.D=D(i);
        %%River incision parameters
        p.Kw =K(j);
        222
        %%
        chi2v=0;
        for k=1:length(DEMv);
            DEM=DEMv{k};
            named=namev{k};
            %%
            p.fileprefix=[named, 'RunK',num2str(p.Kw), 'D', num2str(p.D), 't']; ;
            
            %% Initialize parameter structure.
            % By making p an instance of ttlemset, the user ensures parameter values
            % are set in the right way
            p   = ttlemset(p);
            %% >run the model to a recent amount of time where i have volume of sed
            % computed
            U=DEM;
            U.Z=DEM.Z*0;
            ttlem_out = ttlem(DEM,U,p);
            333
            %% > compute the volume of seds produced - teh total change in elevation
            % within basin and not in lake*area
            out_name=[p.resultsdir named, 'RunK',num2str(p.Kw), 'D', num2str(p.D),'t', num2str(t0), '.mat'];
            load(out_name,'H1');
            h0=H1;
            
            %%
            close all
            dzdt=DEM.Z-h0.Z;
            dzdtb=dzdt*NaN;
            B=Bv{k};
            lk=lkv{k};
            dzdtb(B.Z==1 & lk~=1)=dzdt(B.Z==1 & lk~=1);
            
            %         subplot(1,2,1)
            %         imagesc(DEM),shg
            %         subplot(1,2,2)
            %         imagesc(dzdtb),shg
            %         axis image;
            
            %% > compute teh eroded volume
            Vmod=nansum(dzdtb(:))*pa
            
            
            % >compute the chi^2 between modeled and observed, for the error, use the
            % v/t err. Which is:
            % the
            
            chi2v(k)=(Vmod-Vemp(k)).^2/(Vemp_u(k)^2);
        end
        
        chi2M(i,j)=sum(chi2v);%chi sq sum
    end
end


%

%close all, contour(K,D,log10(chi2M),20),shg
imagesc(Mk(1,:), Md(:,1),log10(chi2M)),shg
xlabel('K [m^{1-2m}yr^{-1}]');
ylabel('D [m^2yr^{-1}]');

%%
%save('chi2Mjanlost_Ath2e5', 'chi2M', 'Md', 'Mk');
% save('chi2Mjanlost_Ath0', 'chi2M', 'Md', 'Mk');

%%
clc
%load 'chi2Mjanlost_Ath2e5'%this is with thresh area
load 'chi2Mjanlost_Ath0'%this is without thresh Area 
%% increase res for visualization
close all
chi2i=(interp2(log10(chi2M),3));
Mdi=(interp2(Md,3));
Mki=(interp2(Mk,3));

%%smooth chi2M
hsize=10;%15
sigma=2;%2
h = fspecial('gaussian', hsize, sigma);
C=imfilter(chi2i, h, 'replicate');

%%plot
imagesc(Mki(1,:), Mdi(:,1), C);
shg

%xlabel('K [L^{1-2m}t^{-1}]');
xlabel('K [m^{1-2m}yr^{-1}]');
ylabel('D [m^2yr^{-1}]');
%set(gca, 'ytick', [0:0.1:0.3], 'xtick', [0:0.002:0.01], 'ydir', 'normal')
hc=colorbar;
title(hc, 'log(\chi^2)')
shg

%% export fig
pathname='/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/figs/';
%fname='all_lakes_rose_pairs';
%fname='all_lakes_pairs';
fname='KD_optim_lost_jan_myr';
name=[pathname, fname, ''];
export_fig(name, '-pdf', '-transparent');

%% %%%%%%% make plot of erosion on topo
%%Plot the dzdt in each of teh nwks
k=2;



%compute the erosion with teh optimal D and K in interp data
[r,c]=find(C==min(C(:)));
ko=Mki(r,c);
do=Mdi(r,c);

%%compute the erosion with teh optimal D and K in otiginal coarse data
% [r,c]=find(chi2M==min(chi2M(:)));
% ko=Mk(r,c);
% do=Md(r,c);



p = ttlemset; %#ok show parameters<NOPTS>

%%run teh model with best k and D


%%Temporal domain
p.TimeSpan=t0;
p.TimeStep=t0;


p.D=do;
%%River incision parameters
p.Kw =ko;

%%Hillslope processes, parameters values
p.m  = 0.5;
p.n  = 1;
p.AreaThresh = 2e5;%0; % channel contributing area threshold, m²

% (un)comment to get the drainage development type of choice
% Fixed or variable drainage netwrok through time
 p.DrainDir='variable';
% p.DrainDir='fixed';


%%Threshold slopes
p.Sc=1; %%21°)
p.Sc_unit='tangent';


%%Output
p.ploteach=inf;
p.saveeach=1;
p
%%Initialize parameter structure.
% By making p an instance of ttlemset, the user ensures parameter values
% are set in the right way
namev={'jan', 'lost'};
named=namev{k};
DEM=DEMv{k};
U=DEM;
U.Z=DEM.Z*0;

p.fileprefix=[named, 'RunK',num2str(p.Kw), 'D', num2str(p.D), 't'];             
p=ttlemset(p);

%





%%> compute the volume of seds produced - teh total change in elevation
% within basin and not in lake*area
ttlem_out = ttlem(DEM,U,p);
out_name=[p.resultsdir named, 'RunK',num2str(p.Kw), 'D', num2str(p.D),'t', num2str(t0), '.mat'];
load(out_name,'H1');
h0=H1;

%
close all
dzdt=DEM.Z-h0.Z;
dzdtb=dzdt*NaN;
B=Bv{k};
lk=lkv{k};

% %%smooth dzdt for removing teh stripe from the dem.
% hsize=10;%15
% sigma=5;%2
% h = fspecial('gaussian', hsize, sigma);
% dzdtb=imfilter(dzdtb, h, 'replicate');

dzdtb(B.Z==1 & lk~=1)=dzdt(B.Z==1 & lk~=1);
vmod=nansum(dzdtb(:)*pa)
testchi=(Vemp(k)-vmod)^2/(Vemp_u(k))^2%vemp-vmod=4.6951e4, 3.29e4 with the interp data - in both the diff is about 0.1 of teh Vemp
min(chi2M(:))
(Vemp(k)-vmod)/Vemp(k)



%% plot DEM+contours using hillshade
% close all, clc
% %imageschs(DEM,[],'colormap',flipud(copper));
% hs=hillshade(DEM);
% imagesc(hs.Z),shg;
% colormap gray
% freezeColors
%%

%hold on
close all;
colormap winter%spring%default%cool;%default; %copper%spring%HSV%jet
colormap (flipud(colormap))%cool;%default; %copper%spring%HSV%jet

dzdtc=dzdtb;
nprc=5;
cthmx=prctile(dzdtb(dzdtb>0),100-nprc);
cthmn=prctile(dzdtb(dzdtb<0),nprc);
dzdtc(dzdtb>cthmx)=cthmx;
dzdtc(dzdtb<cthmn)=cthmn;

%dzdtc_s=dzdtc/range(dzdtc(:))*range(hs(:));
%h=imagesc(-dzdtc_s);%note the minus sigh to illustrate the erosion
h=imagesc(-dzdtc);%note the minus sigh to illustrate the
[r,c]=find(B.Z==1);
np=50;
xlm=[min(c)-np, max(c)+np];
xlim(xlm);
ylm=[min(r)-np, max(r)+np]
ylim(ylm);
shg
%
%I=imadjust(dzdtc);
%I=imadjust(dzdtc,[],[],1);
%h=imagesc(I);
%colorbar(peer,h)
freezeColors
%
set(h, 'AlphaData', ~isnan(dzdtc))
shg
%
hold on
nc=round(range(DEM.Z(:))/10);
contour(DEM.Z, 20, '-k')
set(gca, 'ytick', [], 'xtick', [])
shg
caxis([min(dzdtc(:)), max(dzdtc(:))])
hc=colorbar;
title(hc, 'E[m]')

%%scale bar
xlen=20;
x0=min(xlm)+np/3;
y0=max(ylm)-np/4;
plot([x0,x0+xlen], [y0,y0], '-k', 'linewidth', 2)
text(x0, y0-10, [num2str(xlen*DEM.cellsize), 'm'])

%% explrt fig
pathname='/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/figs/';
%fname='all_lakes_rose_pairs';
%fname='all_lakes_pairs';
fname='lost_erosion';
%fname='jan_erosion';
name=[pathname, fname, ''];
export_fig(name, '-pdf', '-transparent');


%%
close all
imagesc(dzdtb, [-1,1]),shg,colorbar

%%
dzdtc=dzdtc*40;
%dzdtc=dzdtc./range(dzdtc(:))*3;
sh=scatter3(c,r,[],dzdtc(:),'.');
range(dzdtc(:))
colorbar
%
hold on
contour(DEM.Z, 20, '-k')

shg


%% plot the erosion value over the previous plot
close all
hist(dzdtc(:)),shg
%%
close all
imagesc(dzdtc),shg
colorbar


%%plot 
