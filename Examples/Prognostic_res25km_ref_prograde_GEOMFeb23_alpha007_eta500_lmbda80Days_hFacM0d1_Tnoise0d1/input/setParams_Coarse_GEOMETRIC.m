%%%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.
%%% attention: simTime, nIter0, diag_freq_avg, diag_freq_inst, init_extern,
%%% tau_x0, Topographic parameters, tAlpha, beta, alpha gamma, data.pkg, useRBCsalt
%%% ivdc_kappa (0 for kpp or 100 for normal run), MXLDEPTH, Canyon,
%%% tauRelaxT, prograde,  Time step size, bottomDragQuadratic, Domain size in y
%%% GM_isopycK, GM_background_K, surface heat flux, pChkptFreq
function nTimeSteps = setParams (inputpath,codepath,listterm,Nx,Ny,Nr)  
  

  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%  
  
  tau_x0 = 0.05*(1) %%% Wind stress magnitude
    
  Ws = 50*1000*(1)  %%% Slope half-width
 
  tAlpha = 1e-4*(1) %%% Linear thermal expansion coeff
  
  f0 = 1e-4*(1); %%% Coriolis parameter
  
  
  hFacMin = 0.1;
 

%%%%% CHECK diag_fields_inst %%%%%%%
  
  
  %%% Tapering scheme
  GM_taper_scheme     = 'clipping'; % clipping dm95 linear gkw91 ldd97
  %%% Maximum isopycnal slope 
  GM_maxSlope         = 0.1; % Default:0.01  Mak18:0.005  Mak17:0.05  Andrew:0.025  
 
  if strcmp(GM_taper_scheme,'dm95')
  %%% DM95 critical slope
  GM_Scrit            = GM_maxSlope;
  %%% DM95 tapering width
  GM_Sd               = GM_maxSlope/10;
  end
  
  
  % smaller values (e.g., 0.01) will destroy the simulations with continental
  % slopes for most of tapering schemes
  
  %%% Set true if initialization files are set externally
  init_extern = false

  
  %%% Horizontal viscosity  
  viscAh = 0;  
  %%% Grid-dependent biharmonic viscosity
  viscA4Grid = 0.1; 
  %%% non-dimensional Smagorinsky biharmonic viscosity factor
  viscC4smag = 0;
  %%% Vertical viscosity
  viscAr = 0;     % Si22:3e-4  % NeverWorld2: 1e-4
  %%% Vertical temp diffusion
  diffKrT = 0;    % Mak18, Si22  1e-5
   
  
  %%% Random noise amplitude, only works for   init_extern = false
  tNoise = 0.1;   %0.1
  Lx = 800*1000; %%% Domain size in x 
  Ly = 600*1000; %%% Domain size in y  

  
  Prograde = true;   
  GM_background_K   = 0;  
  GM_AdvForm        = false;
  
  
  
  
  GM_useGEOM        = true;
  
  UseMLTrhines       = false;
  
  
  GEOM_SlopeAware   = false;
  
 
  GEOM_alpha        = 0.07; 
%   For UeLrh: 0.05;
%   For GEOMETRIC: 0.07;
  
  UseVarLmbda       = false;
  HRS_Name          = 'ref';

  GEOM_lmbda        = 1.45e-7;  %1.45e-7


 % 11.6 5.79 2.89 for 10 20 40 days
 % 1.93 1.45 1.16 for 60 80 100 days
 % 0.965 0.827 0.723  for 120 140 160 days
 % 0.681 0.643  for 170 180 days
 % 1.65  for 70 days
 
  GEOM_SloAwa_Burger= true;
  GEOM_minval_K     = 0;
  GEOM_maxval_K     = 1.e+4;
  ene_init          = 1.e-2;  %Initial depth-int energy 
  
  ene_const        =  false;
  GEOM_ConstEner   =  0;
  
  ene_kappa         = 500;  %1000
  
  ene_local         = true;
  vert_struc        = false;
  GEOM_DepthTaper   = false;


  
  GEOM_pickup_write_mdsio = true;
  GEOM_pickup_write_mnc   = false;
  GEOM_pickup_read_mdsio  = false;
  GEOM_pickup_read_mnc    = false;
  
  

    
%%
    
    
    
    
    
    
    



    
    
    
    
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true;      
  fignum = 1;
  ieee='b';
  prec='real*8';
  %%% Data format parameters
  realdigits = 8;
  realfmt=['%.',num2str(realdigits),'e'];
  
  %%% Get parameter type definitions
  paramTypes;
  
  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 
  
  %%% Seconds in one hour
  t1min = 60;
  %%% Seconds in one hour
  t1hour = 60*t1min;
  %%% Seconds in one day
  t1day = 24*t1hour;
  %%% Seconds in 1 year
  t1year = 365*t1day;  
  %%% Metres in one kilometre
  m1km = 1000;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  simTime = 120* t1year + 0.5*t1day ; % %%% Simulation time  80* t1year 1* t1year + 
  diag_phase_avg = 60*t1year;    
  
  nIter0 = 0; %%% Initial iteration 

  
  Ls = 50*m1km; %%% Southern sponge width
  Ln = 150*m1km; %%% Northern sponge width
  Lsurf = 10; %%% Surface sponge width Yan
  H = 4000; %%% Domain size in z 
  g = 9.81; %%% Gravity
 
  beta = 0; %4e-11;  %%% Beta parameter      
  

  viscA4 = 0; %%% Biharmonic viscosity
  viscAhGrid = 0; %%% Grid-dependent harmonic viscosity
  diffKhT = 0; %%% Horizontal temp diffusion

  
  
  %%% Topographic parameters
  Zs = 2250;% %%% Vertical slope position
  Hs = 3500;%; %%% Shelf height
  Ys = 200*m1km; %%% Meridional slope position
  
  
  %%% Stratification parameters
  
  %%% Small domain, delta > 0 wind
  use_wind = true;
  restore_south = false;
  mintemp = 0;
  maxtemp = 10;
  difftemp = 5;
  

  %%% Parameters related to periodic wind forcing
  %%% REF: http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-November/008068.html
  %%% REF: http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-November/008064.html
  %%% REF: http://mitgcm.org/public/r2_manual/latest/online_documents/node104.html
  periodicExternalForcing = false;
  externForcingCycle = t1year;
  nForcingPeriods = 50;
  if (~periodicExternalForcing)
    nForcingPeriods = 1;   
  end
  externForcingPeriod = externForcingCycle/nForcingPeriods;  
  
  %%% PARM01
  %%% momentum scheme
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
  %%% viscosity  
  parm01.addParm('viscAr',viscAr,PARM_REAL);
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',true,PARM_BOOL);
  parm01.addParm('viscC4leith',0,PARM_REAL);
  parm01.addParm('viscC4leithD',0,PARM_REAL); 
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL); 
  parm01.addParm('viscC4smag',viscC4smag,PARM_REAL);  
  %%% diffusivity
  parm01.addParm('tempAdvScheme',80,PARM_INT);
  parm01.addParm('saltAdvScheme',80,PARM_INT);
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL);
  parm01.addParm('diffK4T',0,PARM_REAL);  
  parm01.addParm('tempStepping',true,PARM_BOOL);
  parm01.addParm('saltStepping',false,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
  %%% equation of state
  parm01.addParm('eosType','LINEAR',PARM_STR); 
  parm01.addParm('tAlpha',tAlpha,PARM_REAL); 
  parm01.addParm('sBeta',0,PARM_REAL); 
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
  parm01.addParm('bottomDragLinear',0,PARM_REAL);
  parm01.addParm('bottomDragQuadratic',2.5e-3,PARM_REAL); % 2.5e-3 for normal runs, 0 for spindown
  %%% physical parameters
  parm01.addParm('f0',f0,PARM_REAL);
  parm01.addParm('beta',beta,PARM_REAL);
  parm01.addParm('gravity',g,PARM_REAL);
  %%% full Coriolis force parameters
  parm01.addParm('quasiHydrostatic',false,PARM_BOOL);
  parm01.addParm('fPrime',0,PARM_REAL);
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',100,PARM_REAL);
  parm01.addParm('implicitDiffusion',true,PARM_BOOL);
  parm01.addParm('implicitViscosity',true,PARM_BOOL);
  %%% exact volume conservation
  parm01.addParm('exactConserv',true,PARM_BOOL);
  %%% C-V scheme for Coriolis term
  parm01.addParm('useCDscheme',false,PARM_BOOL);
  %%% partial cells for smooth topography
  parm01.addParm('hFacMin',hFacMin,PARM_REAL);  
  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);
  parm01.addParm('debugLevel',1,PARM_INT);
  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL);
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL);


 
  %%% PARM02
  parm02.addParm('useSRCGSolver',true,PARM_BOOL);  
  parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL);
 
  %%% PARM03
  parm03.addParm('alph_AB',1/2,PARM_REAL);
  parm03.addParm('beta_AB',5/12,PARM_REAL);
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
  parm03.addParm('chkptFreq',0*t1year,PARM_REAL);
  parm03.addParm('pChkptFreq',5*t1year,PARM_REAL);
  parm03.addParm('taveFreq',0,PARM_REAL);
  parm03.addParm('dumpFreq',0,PARM_REAL);
  parm03.addParm('monitorFreq',5*t1year,PARM_REAL);
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL);
  if (periodicExternalForcing)
    parm03.addParm('periodicExternalForcing',periodicExternalForcing,PARM_BOOL);
    parm03.addParm('externForcingPeriod',externForcingPeriod,PARM_REAL);
    parm03.addParm('externForcingCycle',externForcingCycle,PARM_REAL);
  end
  
  %%% PARM04
  parm04.addParm('usingCartesianGrid',true,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',false,PARM_BOOL);    
  
    
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GRID SPACING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%    
    

  %%% Zonal grid
  dx = Lx/Nx;  
  xx = (1:Nx)*dx;
  xx = xx-mean(xx);
  
  %%% Uniform meridional grid   
  dy = (Ly/Ny)*ones(1,Ny);  
  yy = cumsum((dy + [0 dy(1:end-1)])/2);
 
  %%% Plotting mesh
  [Y,X] = meshgrid(yy,xx);
  
  %%% Grid spacing increases with depth, but spacings exactly sum to H
  zidx = 1:Nr;
  
  if(Nr == 70)
  gamma = 10;  
  alpha = 10;
  elseif(Nr == 133)
  gamma = 20; 
  alpha = 5; 
  elseif(Nr == 35)
  gamma = 5; 
  alpha = 20;             
  end
            
  dz1 = 2*H/Nr/(alpha+1);
  dz2 = alpha*dz1;
  dz = dz1 + ((dz2-dz1)/2)*(1+tanh((zidx-((Nr+1)/2))/gamma));
  zz = -cumsum((dz+[0 dz(1:end-1)])/2);

  %%% Store grid spacings
  parm04.addParm('delX',dx*ones(1,Nx),PARM_REALS);
  parm04.addParm('delY',dy,PARM_REALS);
  parm04.addParm('delR',dz,PARM_REALS);      
  
  %%% Don't allow partial cell height to fall below min grid spacing
  parm01.addParm('hFacMinDr',min(dz),PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
   

  %%% tanh-shaped slope  
  h = - Zs - (Hs/2)*tanh((Y-Ys)/Ws);
%==========================================================================
  %%% Add Canyons by Yan
%   Ys_curve = Ys + (2)*(Ws/2)*cos(2*pi*5/Lx*xx);
%   for i = 1:Nx
%       h(i, :) = - Zs - (Hs/2)*tanh((Y(i, :)-Ys_curve(i))/Ws);
%   end

%   %%% Add troughs
%   np = 5;
%   ycs = (Ys - Ws/4) + (Ws/2)*cos(2*pi*np/Lx*xx);
%   yns = (Ys + Ws/4) + (Ws/2)*cos(2*pi*np/Lx*xx);
%   for i = 1:Nx
%       ind = and(yy > ycs(i), yy < yns(i));
%       h(i, ind) = -H;
%   end
%==========================================================================  
  h(:,1) = 0;   
  h(:,end) = 0;  
  

  
  %%% Plot the bathymetry
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(Y(1,:),h(1,:));
    title('Model bathymetry');
  end
  disp(['shelf depth ' num2str(h(1, 2))]) 
  
  %%% Save as a parameter
  writeDataset(h,fullfile(inputpath,'bathyFile.bin'),ieee,prec);
  parm05.addParm('bathyFile','bathyFile.bin',PARM_STR); 
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% NORTHERN TEMPERATURE/SALINITY PROFILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Exponential background temperature  
  Hscale = 1000;  
  tNorth = mintemp + (maxtemp-mintemp)*(exp(zz/Hscale)-exp(-H/Hscale))/(1-exp(-H/Hscale));
  if (restore_south)
    if (Prograde)
      tSouth = maxtemp + difftemp + zz*N2South/g/tAlpha;
      tSouth(zz<-(H-Hs)) = mintemp + (maxtemp + difftemp - (H-Hs)*N2South/g/tAlpha)*(zz(zz<-(H-Hs))+H)/Hs;
    else
      tSouth = mintemp + difftemp + (zz+(H-Hs))*N2South/g/tAlpha;
      tSouth(zz<-(H-Hs)) = mintemp + difftemp*(zz(zz<-(H-Hs))+H)/Hs;
    end
  else
    tSouth = tNorth;
  end
  
  %%% Constant background salinity (currently unused)
  sRef = 1;    
  sNorth = sRef*ones(1,Nr);
  sSouth = 0*ones(1,Nr);
  
  %%% Plot the northern temperature
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(tNorth,zz);
    xlabel('\theta_n_o_r_t_h');
    ylabel('z','Rotation',0);
    title('Northern relaxation temperature');
  end
  
  %%% Plot the southern temperature
  if (restore_south)
    if (showplots)
      figure(fignum);
      fignum = fignum + 1;
      clf;
      plot(tSouth,zz);
      xlabel('T_s_o_u_t_h');
      ylabel('z','Rotation',0);
      title('Southern relaxation temperature');
    end
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Weights for linear interpolation
  wtNorth = repmat(Y/Ly,[1 1 Nr]);
  wtSouth = 1 - wtNorth;
      
  %%% Linear variation between southern and northern reference temperatures
  hydroThNorth = repmat(reshape(tNorth,[1,1,Nr]),[Nx Ny]);
  hydroThSouth = repmat(reshape(tSouth,[1,1,Nr]),[Nx Ny]); 
  hydroTh = wtNorth.*hydroThNorth + wtSouth.*hydroThSouth;
  TT = squeeze(hydroTh(1,:,:));
  
  %%% Plot the initial condition
  if (showplots)
  
    
     figure(fignum);
    fignum = fignum + 1;
    clf;
    [ZZ,YY]=meshgrid(zz,yy);
    hydroTh_plot = squeeze(hydroTh(1,:,:));
    HH = repmat(reshape(h(1,:),[Ny,1]),[1 Nr]);
    hydroTh_plot(ZZ<HH) = NaN;
%     contourf(YY,ZZ,hydroTh_plot,30);
    pcolor(YY,ZZ,hydroTh_plot);
    colorbar;
    xlabel('y');
    ylabel('z');
    title('Initial temperature section');
    
    
  end
  
  %%% Linear initial salinity
  hydroSaNorth = repmat(reshape(sNorth,[1,1,Nr]),[Nx Ny]);
  hydroSaSouth = repmat(reshape(sSouth,[1,1,Nr]),[Nx Ny]);   
  hydroSa = wtNorth.*hydroSaNorth + wtSouth.*hydroSaSouth; 
    
  %%% Add some random noise
  hydroTh = hydroTh + tNoise*(2*rand(Nx,Ny,Nr)-1);  
    
  %%% Write to data files
  %%% N.B. These will not be used because the files will be overwritten by
  %%% those in the DEFAULTS folder
  if (~init_extern)
    writeDataset(hydroTh,fullfile(inputpath,'hydrogThetaFile.bin'),ieee,prec); 
    parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);  
%     writeDataset(hydroSa,fullfile(inputpath,'hydrogSaltFile.bin'),ieee,prec); 
%     parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR);  
  else
    parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);
%     parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR);  
    parm05.addParm('uVelInitFile','uVelInitFile.bin',PARM_STR);  
    parm05.addParm('vVelInitFile','vVelInitFile.bin',PARM_STR);  
    parm05.addParm('pSurfInitFile','pSurfInitFile.bin',PARM_STR);  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFORMATION RADIUS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %%% Temperature stratification   
  DZ = repmat(dz,[Ny 1]);
  Tz = zeros(size(TT));  
  Tz(:,1) = (TT(:,1)-TT(:,2)) ./ (0.5*(DZ(:,1)+DZ(:,2)));
  Tz(:,2:Nr-1) = (TT(:,1:Nr-2)-TT(:,3:Nr)) ./ (DZ(:,2:Nr-1) + 0.5*(DZ(:,1:Nr-2)+DZ(:,3:Nr)));
  Tz(:,Nr) = (TT(:,Nr-1)-TT(:,Nr)) ./ (0.5*(DZ(:,Nr-1)+DZ(:,Nr)));  
  N2 = tAlpha*g*Tz;

  %%% Calculate internal wave speed and first Rossby radius of deformation
  N = sqrt(N2);
  Cig = zeros(size(yy));
  for j=1:Ny    
    for k=1:Nr
      if (zz(k) > h(1,j))        
        Cig(j) = Cig(j) + N(j,k)*dz(k);
      end
    end
  end
  Rd = Cig./(pi*abs(f0+beta*Y(1,:)));

  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    semilogx(N2,zz/1000);
    xlabel('N^2 (km)');
    ylabel('z (km)','Rotation',0);
    title('Buoyancy frequency');
  end
  
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy/1000,Rd/1000);
    xlabel('y (km)');
    ylabel('R_d (km)','Rotation',0);
    title('First baroclinic Rossby deformation radius');
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATE TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  
  %%% These estimates are in no way complete, but they give at least some
  %%% idea of the time step needed to keep things stable. In complicated 
  %%% simulations, preliminary tests may be required to estimate the
  %%% parameters used to calculate these time steps.        
  
  %%% Gravity wave CFL

  %%% Upper bound for absolute horizontal fluid velocity (m/s)
  %%% At the moment this is just an estimate
  Umax = 2.0;  
  %%% Max gravity wave speed 
  cmax = max(Cig);
  %%% Max gravity wave speed using total ocean depth
  cgmax = Umax + cmax;
  %%% Advective CFL
  deltaT_adv = min([0.5*dx/cmax,0.5*dy/cmax]);
  %%% Gravity wave CFL
  deltaT_gw = min([0.5*dx/Umax,0.5*dy/Umax]);
  %%% CFL time step based on full gravity wave speed
  deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);
    
  %%% Other stability conditions
  
  %%% Inertial CFL time step (Sf0<=0.5)
  deltaT_itl = 0.5/abs(f0);
  %%% Time step constraint based on horizontal diffusion 
  deltaT_Ah = 0.5*min([dx dy])^2/(4*viscAh);    
  %%% Time step constraint based on vertical diffusion
  deltaT_Ar = 0.5*min(dz)^2 / (4*viscAr);  
  %%% Time step constraint based on biharmonic viscosity 
  deltaT_A4 = 0.5*min([dx dy])^4/(32*viscA4);
  %%% Time step constraint based on horizontal diffusion of temp 
  deltaT_KhT = 0.4*min([dx dy])^2/(4*diffKhT);    
  %%% Time step constraint based on vertical diffusion of temp 
  deltaT_KrT = 0.4*min(dz)^2 / (4*diffKrT);
  
  %%% Time step size  
  deltaT = min([deltaT_fgw deltaT_gw deltaT_adv deltaT_itl deltaT_Ah deltaT_Ar deltaT_KhT deltaT_KrT deltaT_A4]);
  deltaT = round(deltaT)
  
  
  
%     deltaT = 1000;
  

  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT;
  
  %%% Write end time time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL); 
    
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% ZONAL WIND %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
  
    
 
  %%% Set up a wind stress profile    
  if (use_wind)
    Lwind = 400*m1km;
    tau_x = tau_x0.*sin(pi*yy/Lwind).^2;
    tau_x(yy>Lwind) = 0;
    if (Prograde)
       tau_x = abs(tau_x);
    else
       tau_x = -abs(tau_x);
     end
  else
    tau_x = 0*yy;
  end

  
  %%% Create wind stress matrix  
  tau_mat = zeros(size(Y)); %%% Wind stress matrix  
  for j=1:Ny   
    tau_mat(:,j) = tau_x(j);          
  end 
  
  %%% Plot the wind stress 
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy,squeeze(tau_mat(1,:)));
    xlabel('y');
    ylabel('\tau_w','Rotation',0);
    title('Wind stress profile');
  end
  
  %%% Save as a parameter  
  writeDataset(tau_mat,fullfile(inputpath,'zonalWindFile.bin'),ieee,prec); 
  parm05.addParm('zonalWindFile','zonalWindFile.bin',PARM_STR);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SURFACE HEAT/SALT FLUXES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%%% Heat flux profile
  heat_profile = 0*yy;
  %%% Done
  heat_flux = zeros(Nx,Ny,nForcingPeriods);
  for n=1:nForcingPeriods
    for j=1:Ny      
      heat_flux(:,j,n) = heat_profile(j);             
    end         
  end

  %%% Plot the surface heat flux
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy,squeeze(heat_flux(1,:,:)));
    xlabel('y');
    ylabel('heat flux');
    title('Surface heat flux');
  end  
  
  %%% Save as parameters
  writeDataset(heat_flux,fullfile(inputpath,'surfQfile.bin'),ieee,prec);
  parm05.addParm('surfQfile','surfQfile.bin',PARM_STR);  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
  
 
  
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% RBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  rbcs_parm01 = parmlist;
  rbcs_parm02 = parmlist;
  RBCS_PARM = {rbcs_parm01,rbcs_parm02};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  useRBCtemp = true;
  useRBCsalt = false;
  useRBCuVel = false;
  useRBCvVel = false;
    
if (Prograde)
  tauRelaxT = 14*t1day; % 7d for normal runs, 14d for prograde run
else
  tauRelaxT = 7*t1day; % 7d for normal runs, 14d for prograde run
end
  tauRelaxS = t1day;
  tauRelaxU = t1day;
  tauRelaxV = t1day;
  rbcs_parm01.addParm('useRBCtemp',useRBCtemp,PARM_BOOL);
  rbcs_parm01.addParm('useRBCsalt',useRBCsalt,PARM_BOOL);
  rbcs_parm01.addParm('useRBCuVel',useRBCuVel,PARM_BOOL);
  rbcs_parm01.addParm('useRBCvVel',useRBCvVel,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxT',tauRelaxT,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxS',tauRelaxS,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxU',tauRelaxU,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxV',tauRelaxV,PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION TEMPERATURE/SALINITY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Relaxation T/S matrices
  temp_relax = ones(Nx,Ny,Nr);  
%   salt_relax = ones(Nx,Ny,Nr);  
  Ny_south = round(Ny/2);
  Ny_north = Ny - Ny_south;
  temp_relax(:,1:Ny_south,:) = repmat(reshape(tSouth,[1 1 Nr]),[Nx Ny_south 1]);  
  temp_relax(:,Ny_south+1:Ny,:) = repmat(reshape(tNorth,[1 1 Nr]),[Nx Ny_north 1]);  
%   salt_relax(:,1:Ny_south,:) = repmat(reshape(sSouth,[1 1 Nr]),[Nx Ny_south 1]);  
%   salt_relax(:,Ny_south+1:Ny,:) = repmat(reshape(sNorth,[1 1 Nr]),[Nx Ny_north 1]);
  
  %%% Yan: relax SST linearly from 15(onshore) to 10(offshore) for prograde
  %%% run
  if (Prograde)
      
  indy = find(yy <= 4.5e5, 1, 'last');
  Tsurf = linspace(15, 10, indy+1);
%   for j = 1:indy+1
   for j = 1:indy     
      temp_relax(:, j, 1:10) = Tsurf(j);
  end
  
  end
  
  
  %%% Plot the relaxation temperature
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    [ZZ,YY]=meshgrid(zz,yy);
%     contourf(YY,ZZ,squeeze(temp_relax(1,:,:)),30);
    pcolor(YY,ZZ,squeeze(temp_relax(1,:,:)));
    colorbar;
    xlabel('y');
    ylabel('z');
    title('Relaxation temperature');
  end
  
  %%% Plot the relaxation salinity
%   if (showplots)
%     figure(fignum);
%     fignum = fignum + 1;
%     clf;
%     [ZZ,YY]=meshgrid(zz,yy);
%     contourf(YY,ZZ,squeeze(salt_relax(1,:,:)),30);
%     colorbar;
%     xlabel('y');
%     ylabel('z');
%     title('Relaxation salinity');
%   end
  
  %%% Save as parameters
  writeDataset(temp_relax,fullfile(inputpath,'sponge_temp.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxTFile','sponge_temp.bin',PARM_STR);  
  
%   writeDataset(salt_relax,fullfile(inputpath,'sponge_salt.bin'),ieee,prec); 
%   rbcs_parm01.addParm('relaxSFile','sponge_salt.bin',PARM_STR);  
  
  
  %%%%%%%%%%%%%%%%%%%%%  
  %%%%% RBCS MASK %%%%%
  %%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Mask is zero everywhere by default, i.e. no relaxation
  msk_temp = zeros(Nx,Ny,Nr);  
%   msk_salt = zeros(Nx,Ny,Nr);  
  
  %%% Set relaxation timescales 
  northRelaxFac = 1;  
  southRelaxFac = 1;  
  surfaceRelaxFac = 1; 
  %%% If a sponge BC is required, gradually relax in the gridpoints 
  %%% approaching the wall (no relaxation at the wall)    
  for j=1:Ny
    for k=1:Nr              
      if ( ( and( yy(j)>450e3,yy(j)<500e3)   ) )
        msk_temp(:,j,k) = (1/northRelaxFac) * (yy(j)/1e3-450)/50;      % prograde/spindown comment to remove north sponge    
      elseif(yy(j)>500e3)        
        msk_temp(:,j,k) = 1;      
      end
    end
  end
  
if (Prograde)
  for j = 1:indy+1
      for k = 1:Nr
           if zz(k) > 0-Lsurf % Yan: for prograde run
              msk_temp(:,j,k) = (1/surfaceRelaxFac) * (zz(k)-(-Lsurf)) / Lsurf;
           end
      end
  end
end
  %%% Plot the temperature relaxation timescale
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    [ZZ,YY]=meshgrid(zz,yy);
%     contourf(YY,ZZ,squeeze(msk_temp(1,:,:)),30);
    pcolor(YY,ZZ,squeeze(msk_temp(1,:,:)))
    colorbar;
    xlabel('y');
    ylabel('z');
    title('Temperature relaxation fraction');
  end  
  
%   %%% Plot the salinity relaxation timescale
%   if (showplots)
%     figure(fignum);
%     fignum = fignum + 1;
%     clf;
%     [ZZ,YY]=meshgrid(zz,yy);
%     contourf(YY,ZZ,squeeze(msk_salt(1,:,:)),30);
%     colorbar;
%     xlabel('y');
%     ylabel('z');
%     title('Salinity relaxation fraction');
%   end  
   
  %%% Save as an input parameters
  writeDataset(msk_temp,fullfile(inputpath,'rbcs_temp_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_temp_mask.bin',PARM_STR);
%   writeDataset(msk_salt,fullfile(inputpath,'rbcs_salt_mask.bin'),ieee,prec); 
%   rbcs_parm01.addParm('relaxMaskFile(2)','rbcs_salt_mask.bin',PARM_STR);  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.rbcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
  %%% Creates the 'data.rbcs' file
  write_data_rbcs(inputpath,RBCS_PARM,listterm,realfmt);
  
  
    
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% OBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  obcs_parm01 = parmlist;
  obcs_parm02 = parmlist;
  OBCS_PARM = {obcs_parm01,obcs_parm02};  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFINE OPEN BOUNDARY TYPES (OBCS_PARM01) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Enables an Orlanski radiation condition at the northern boundary
  useOrlanskiNorth = true;
  OB_Jnorth = Ny*ones(1,Nx);     
  obcs_parm01.addParm('useOrlanskiNorth',useOrlanskiNorth,PARM_BOOL);
  obcs_parm01.addParm('OB_Jnorth',OB_Jnorth,PARM_INTS); 
  
  %%% Enforces mass conservation across the northern boundary by adding a
  %%% barotropic inflow/outflow
  useOBCSbalance = true;
  OBCS_balanceFacN = -1; 
  OBCS_balanceFacE = 0;
  OBCS_balanceFacS = 0;
  OBCS_balanceFacW = 0;
  obcs_parm01.addParm('useOBCSbalance',useOBCSbalance,PARM_BOOL);  
  obcs_parm01.addParm('OBCS_balanceFacN',OBCS_balanceFacN,PARM_REAL);  
  obcs_parm01.addParm('OBCS_balanceFacE',OBCS_balanceFacE,PARM_REAL);  
  obcs_parm01.addParm('OBCS_balanceFacS',OBCS_balanceFacS,PARM_REAL);  
  obcs_parm01.addParm('OBCS_balanceFacW',OBCS_balanceFacW,PARM_REAL);  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ORLANSKI OPTIONS (OBCS_PARM02) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Velocity averaging time scale - must be larger than deltaT.
  %%% The Orlanski radiation condition computes the characteristic velocity
  %%% at the boundary by averaging the spatial derivative normal to the 
  %%% boundary divided by the time step over this period.
  %%% At the moment we're using the magic engineering factor of 3.
  cvelTimeScale = 3*deltaT;
  %%% Max dimensionless CFL for Adams-Basthforth 2nd-order method
  CMAX = 0.45; 
  
  obcs_parm02.addParm('cvelTimeScale',cvelTimeScale,PARM_REAL);
  obcs_parm02.addParm('CMAX',CMAX,PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.obcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data.obcs' file
  write_data_obcs(inputpath,OBCS_PARM,listterm,realfmt);
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%% LAYERS %%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  layers_parm01 = parmlist;
  LAYERS_PARM = {layers_parm01};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Define parameters for layers package %%%
  
  %%% Number of fields for which to calculate layer fluxes
  layers_maxNum = 1;
  
  %%% Specify potential temperature
  layers_name = 'TH';  
  
  %%% Potential temperature bounds for layers  
   %%% Potential temperature bounds for layers  
  if (Prograde)
%     layers_bounds = mintemp:0.5:maxtemp+difftemp;
    layers_bounds = [mintemp flip(tNorth) maxtemp:0.1:maxtemp+difftemp];
   else
%     layers_bounds = mintemp:(maxtemp-mintemp)/100:maxtemp;
    layers_bounds = [mintemp flip(tNorth) maxtemp];
  end
  
  %%% Reference level for calculation of potential density
  layers_krho = 1;    
  
  %%% If set true, the GM bolus velocity is added to the calculation
  layers_bolus = false;  
   
  %%% Layers
  layers_parm01.addParm('layers_bounds',layers_bounds,PARM_REALS); 
  layers_parm01.addParm('layers_krho',layers_krho,PARM_INT); 
  layers_parm01.addParm('layers_name',layers_name,PARM_STR); 
  layers_parm01.addParm('layers_bolus',layers_bolus,PARM_BOOL); 

  %%z% Create the data.layers file
  write_data_layers(inputpath,LAYERS_PARM,listterm,realfmt);
  
  %%% Create the LAYERS_SIZE.h file
  createLAYERSSIZEh(codepath,length(layers_bounds)-1,layers_maxNum); 
  
  
  
  
    
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%% GMREDI %%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GMREDI SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  gmredi_parm01 = parmlist;
  GMREDI_PARM = {gmredi_parm01};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GMREDI PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Define parameters for gmredi package %%%
  %%%http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/010377.html%%%

  %%% Isopycnal diffusivity
  GM_isopycK          = 0;







  %%% Add parameters
  gmredi_parm01.addParm('GM_isopycK',GM_isopycK,PARM_REAL);
  gmredi_parm01.addParm('GM_background_K',GM_background_K,PARM_REAL);
  gmredi_parm01.addParm('GM_maxSlope',GM_maxSlope,PARM_REAL);
  gmredi_parm01.addParm('GM_taper_scheme',GM_taper_scheme,PARM_STR);
  gmredi_parm01.addParm('GM_AdvForm',GM_AdvForm,PARM_BOOL);


  
  
  
  
  
  if strcmp(GM_taper_scheme,'dm95')
  gmredi_parm01.addParm('GM_Scrit',GM_Scrit,PARM_REAL);
  gmredi_parm01.addParm('GM_Sd',GM_Sd,PARM_REAL);
  end
  
  
  
  %%% Add GEOMETRIC parameters 


  gmredi_parm01.addParm('vert_struc',vert_struc,PARM_BOOL);
  
%  gmredi_parm01.addParm('GM_SlopeAware',GEOM_SlopeAware,PARM_BOOL);
  gmredi_parm01.addParm('GM_useGEOM',GM_useGEOM,PARM_BOOL);
  
  

  gmredi_parm01.addParm('GEOM_lmbda',GEOM_lmbda,PARM_REAL);
  gmredi_parm01.addParm('GEOM_alpha',GEOM_alpha,PARM_REAL);
  gmredi_parm01.addParm('GEOM_minval_K',GEOM_minval_K,PARM_REAL);
  gmredi_parm01.addParm('GEOM_maxval_K',GEOM_maxval_K,PARM_REAL);
  gmredi_parm01.addParm('ene_local',ene_local,PARM_BOOL);
  gmredi_parm01.addParm('ene_init',ene_init,PARM_REAL);
  
  gmredi_parm01.addParm(' GEOM_ConstEner', GEOM_ConstEner,PARM_REAL);
  
 
  gmredi_parm01.addParm('ene_kappa',ene_kappa,PARM_REAL);
  gmredi_parm01.addParm('GEOM_DepthTaper',GEOM_DepthTaper,PARM_BOOL);
  gmredi_parm01.addParm('GEOM_SlopeAware',GEOM_SlopeAware,PARM_BOOL);
  gmredi_parm01.addParm('GEOM_SloAwa_Burger',GEOM_SloAwa_Burger,PARM_BOOL);
  gmredi_parm01.addParm('UseVarLmbda',UseVarLmbda,PARM_BOOL);
  gmredi_parm01.addParm('ene_const',ene_const,PARM_BOOL);
  gmredi_parm01.addParm('UseMLTrhines',UseMLTrhines,PARM_BOOL);

  
   gmredi_parm01.addParm('GEOM_pickup_write_mdsio',GEOM_pickup_write_mdsio,PARM_BOOL);
   gmredi_parm01.addParm('GEOM_pickup_write_mnc',GEOM_pickup_write_mnc,PARM_BOOL);
   gmredi_parm01.addParm('GEOM_pickup_read_mdsio',GEOM_pickup_read_mdsio,PARM_BOOL);
   gmredi_parm01.addParm('GEOM_pickup_read_mnc',GEOM_pickup_read_mnc,PARM_BOOL);

  
 
 if(UseVarLmbda)
     
   
   eval( ['load G:\TracerScripts\2D_prograde\KgmThetaUvel_Prog_' HRS_Name '_5Ensemble.mat DissipationRate YY2D_ref']);

  GEOM_Varlmbda = coarse25(YY2D_ref(:,1),DissipationRate',yy/1000);
  
  mm=find(isnan(GEOM_Varlmbda));
  GEOM_Varlmbda(mm(1))=GEOM_Varlmbda(mm(1)+1);
  GEOM_Varlmbda(mm(2:end))=GEOM_Varlmbda(mm(2)-1);
  GEOM_Varlmbda = repmat(reshape(GEOM_Varlmbda,[1 Ny]),[Nx 1]);


   writeDataset(GEOM_Varlmbda,fullfile(inputpath,'LmbdaField.bin'),ieee,prec);    

 
      
      
    figure;
    clf;
    plot(yy,mean(GEOM_Varlmbda))
    xlabel('y');
    ylabel('Diagnosed disspation rate');  
 end

%%z% Create the data.gmredi file
  write_data_gmredi(inputpath,GMREDI_PARM,listterm,realfmt);
  
 
 

 
 
 
 
 
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
%   




  diag_fields_avg = ...
  {...
    'UVEL','VVEL','WVEL',... %%% Velocities
    'UVELSQ','VVELSQ','WVELSQ',... %%% For Kinetic Energy
    'UV_VEL_Z', 'UV_VEL_C','WU_VEL','WV_VEL'... %%% Momentum fluxes
    'UVELTH','VVELTH','WVELTH', ... %%% Temperature fluxes  
    'THETA','THETASQ', ... %%% Temperature
     'PHIHYD', ... %%% Pressure
     'MXLDEPTH', ...
     'UUU_VEL', 'UVV_VEL', 'UWW_VEL', ...
     'VUU_VEL', 'VVV_VEL', 'VWW_VEL', ...
     'WUU_VEL', 'WVV_VEL', 'WWW_VEL', ...
     'UVELPHI', 'VVELPHI', 'WVELPHI', ...
     'momKE', 'momHDiv', 'momVort3', ...          %%% momentum diagnostics
     'USidDrag', 'VSidDrag', ...
     'Um_Diss', 'Vm_Diss', 'Um_Advec', 'Vm_Advec', ...
     'Um_Cori', 'Vm_Cori', 'Um_Ext', 'Vm_Ext', ...
     'Um_AdvZ3', 'Vm_AdvZ3', 'Um_AdvRe', 'Vm_AdvRe', ...
     'VISCA4D','VISCA4Z', ...
     'Strain','Tension','Stretch', ...
     'botTauX','botTauY', ...
     'Um_dPhiX','TOTUTEND', ...
     'GEOMkap0','GEOMeE','GEOMEgen','GEOMEdis','GEOMEadv',...
     'GEOMElap',...
     'GM_GEOMK',...
     'SBurger', ...
%    'GEOMEwav','GEOM_c1','GEOMcros','depthC','d_taper',
%    'GEOMstru',
%    'SDelta',...
  };

%     'UVELSLT','VVELSLT','WVELSLT', ... %%% Salt fluxes
%     'SALT','SALTSQ', ... %%% Salinity
%      'ADVy_SLT', 'ADVr_SLT', 'DFyE_SLT', 'DFrE_SLT', 'DFrI_SLT', ...


%      'ADVx_Um', 'ADVy_Um', 'ADVrE_Um', 'ADVx_Vm', ...
%      'ADVy_Vm', 'ADVrE_Vm', 'VISCx_Um', 'VISCy_Um', ...
%      'VISrE_Um', 'VISrI_Um', 'VISCx_Vm', 'VISCy_Vm', ...
%      'VISrE_Vm', 'VISrI_Vm', ...
  
  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = 20*t1year;
%  diag_freq_avg = 5*t1year;

     
  diag_parm01.addParm('diag_mnc',false,PARM_BOOL);  
  for n=1:numdiags_avg    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);  
    
  end

  diag_fields_inst = ...
 {...
   'UVEL','VVEL','WVEL','THETA', ...
   'GEOMkap0','GEOMeE','GEOMEgen','GEOMEdis','GEOMEadv',...
   'GEOMElap',...
   'GM_GEOMK',...
   'SBurger', ...
%    'GEOMEwav','GEOM_c1','GEOMcros','depthC','d_taper',
%    'GEOMstru',
%    'SDelta',...
  };



  numdiags_inst = length(diag_fields_inst);  
  diag_freq_inst = 1*t1year;
%   diag_freq_inst = 1;
%   diag_freq_inst = 1*t1day;
  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    ndiags = ndiags + 1;
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);  
  end
  
  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  createDIAGSIZEh(codepath,ndiags,max(Nr,length(layers_bounds)-1));
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates a matlab file defining all input parameters
  write_matlab_params(inputpath,[PARM RBCS_PARM OBCS_PARM DIAG_MATLAB_PARM LAYERS_PARM],realfmt);
 
  
end
