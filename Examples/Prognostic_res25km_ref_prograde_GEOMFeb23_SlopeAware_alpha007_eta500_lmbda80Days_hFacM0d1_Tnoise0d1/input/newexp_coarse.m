%%%
%%% newexp.m
%%%
%%% Convenience script to set up folders and options for a new MITgcm run
%%% in a way that allows:
%%% - automatic generation of basic folder and file structure
%%% - automatic calculation of inter-dependent input parameters
%%% - compatibility between required numbers of processors
%%%
%%% Sets up a new experiment folder with subdirectories for the build, the 
%%% code options, the input files, and a results folder. Creates a 'build'
%%% script in the build folder and a 'run' script in the results folder.
%%% Generates the SIZE.h file in the code folder based on parameters
%%% specified here, and copies other code files from the DEFAULTS/code 
%%% folder. Generates all 'eedata' and some 'data' parameters in the input 
%%% folder based on parameters specified here and code in create_data 
%%% function. Other parameters are copied from the DEFAULTS/input folder.
%%%
%%% NOTE: 'expname' MUST NOT be set to 'DEFAULTS'
%%% check setParams


clc
close all
clear all


%runsdir = '/home/hweiaf/huaiyu/MITgcm';
runsdir = 'H:\Data\MITgcm_2DSlope\prograde';
OLy = 3; %%% no. of overlapping y-gridpoints per tile 
 
 

% Prognostic_res25km_ref_Coarse_prograde_MakGEOM_SlopeAware_alpha005_lmbda100Days_cap10000
exp_name = 'Prognostic_res25km_ref_prograde_GEOMFeb23_SlopeAware_alpha007_eta500_lmbda80Days_hFacM0d1_Tnoise0d1'
setParamsNAME = 'setParams_Coarse_GEOMETRIC';
ResNodes   = '25km1Nodes600km';  

% ref 32Wind DoubleWind HalfSlope 0d66Slope 32Slope HalfExpansion DoubleExpansion 0d75Coriolis 1d25Coriolis

% exp_name = 'SM_prod_res25km_ref_prograde_GM300_hFacM0d1_Tnoise0d1'
% setParamsNAME = 'setParams_Coarse';
% ResNodes   = '25km1Nodes600km';  



% exp_name = 'Prognostic_res25km_ref_prograde_Visb_L50km_hFacM0d1_Nr35'
% setParamsNAME = 'setParams_Coarse_Visbeck';
% ResNodes   = '25km1Nodes600kmNr35';  
%%%% NEED to modify ./DEFAULTS/code/GMREDI_OPTIONS

% exp_name = 'SM_prod_res5km_FlatBottomBeta_retrograde_GM0_hFacM0d1'
% setParamsNAME = 'setParams_Coarse_FlatBottom';
% ResNodes   = '4km5Nodes600km';  

% exp_name = 'SM_prod_res4km_2Expansion0d3Coriolis_retrograde_GM0'
% setParamsNAME = 'setParams_Coarse_retrograde';
% ResNodes   = '4km4Nodes';  

addpath  G:\MITgcm\newexp_utils
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% USER-SPECIFIED OPTIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  
  %%% Local directory in which to create experiments

  %%% Experiment subdirectories 
  builddir = 'build';
  codedir = 'code';
  inputdir = 'input';
  resultsdir = 'results';

  %%% List terminator character for parameter files - may be '/' or '&'
  %%% depending on operating system
  listterm = '&';

  %%% Line feed character - important for .sh shell script files
  %%% On unix this probably need to be '\n', in windows '\r\n'
  lf = '\n';    
  
  %Typically nTx=nSx and nTy=nSy,
  %%% so that each process handles nSx*nSy tiles, each with its own thread.
  
  nSx = 1; %%% no. of tiles per processor in x-direction
  nSy = 1; %%% no. of tiles per processor in y-direction
  nTx = 1; %%% no. of threads per processor in x-direction
  nTy = 1; %%% no. of threads per processor in y-direction
  OLx = 3; %%% no. of overlapping x-gridpoints per tile


  
  %%% These parameters are most likely to vary between experiments
  %%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  
  

  
switch(ResNodes)
    
    %512=2*2*2*2*2*2*2*2*2
    %800=2*2*2*2*2*5*5
    %14*40= 2*2*2*2*5*7
    
   case('1km1Nodes')
% %  %%% Set-up for hpc3, 800x512x70, 40 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 160; %%% no. of x-gridpoints per tile
  sNy = 64 ; %%% no. of y-gridpoints per tile
  nPx = 5; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster  
    
  case('1km5Nodes')
% %  %%% Set-up for hpc3, 800x512x70, 200 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 32; %%% no. of x-gridpoints per tile
  sNy = 64 ; %%% no. of y-gridpoints per tile
  nPx = 25; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
    
  
    case('1km5Nodes800km')
% %  %%% Set-up for hpc3, 800x512x70, 200 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 80; %%% no. of x-gridpoints per tile
  sNy = 40; %%% no. of y-gridpoints per tile
  nPx = 10; %%% no. of processors in x-direction
  nPy = 20; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
    
  
    
    case('1km8Nodes')
% %  %%% Set-up for hpc3, 800x512x70, 320 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 40; %%% no. of x-gridpoints per tile
  sNy = 32 ; %%% no. of y-gridpoints per tile
  nPx = 20; %%% no. of processors in x-direction
  nPy = 16; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster


  
  
    case('1km10Nodes')
%   %%% Set-up for hpc3, 800x512x70, 400 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 32; %%% no. of x-gridpoints per tile
  sNy = 32 ; %%% no. of y-gridpoints per tile
  nPx = 25; %%% no. of processors in x-direction
  nPy = 16; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster

  
   case('2km4Nodes')
% %  %%% Set-up for hpc3, 400x256x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 40; %%% no. of x-gridpoints per tile
  sNy = 16 ; %%% no. of y-gridpoints per tile
  nPx = 10; %%% no. of processors in x-direction
  nPy = 16; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
 case('2km5Nodes')
%   %%% Set-up for hpc3, 400x256x70, 200 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 16; %%% no. of x-gridpoints per tile
  sNy = 32 ; %%% no. of y-gridpoints per tile
  nPx = 25; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster
  
   case('2km8Nodes')
% %  %%% Set-up for hpc3, 400x256x70, 320 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 20; %%% no. of x-gridpoints per tile
  sNy = 16 ; %%% no. of y-gridpoints per tile
  nPx = 20; %%% no. of processors in x-direction
  nPy = 16; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster

  
     case('2km5NodesLy600')
% %  %%% Set-up for hpc3, 400x300x70, 200 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 20; %%% no. of x-gridpoints per tile
  sNy = 30 ; %%% no. of y-gridpoints per tile
  nPx = 20; %%% no. of processors in x-direction
  nPy = 10; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
     case('4km4Nodes')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 20; %%% no. of x-gridpoints per tile
  sNy = 8 ; %%% no. of y-gridpoints per tile
  nPx = 10; %%% no. of processors in x-direction
  nPy = 16; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
    
    case('8km1Nodes')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 20; %%% no. of x-gridpoints per tile
  sNy = 8 ; %%% no. of y-gridpoints per tile
  nPx = 5; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
 
  case('16km1Nodes')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 10; %%% no. of x-gridpoints per tile
  sNy = 4 ; %%% no. of y-gridpoints per tile
  nPx = 5; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
    case('32km1Nodes')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 5; %%% no. of x-gridpoints per tile
  sNy = 2 ; %%% no. of y-gridpoints per tile
  nPx = 5; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
  
  case('4km5Nodes600km')
% %  %%% Set-up for hpc3, 200x150x70, 200 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 10; %%% no. of x-gridpoints per tile
  sNy = 15 ; %%% no. of y-gridpoints per tile
  nPx = 20; %%% no. of processors in x-direction
  nPy = 10; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  case('4km4Nodes')
% %  %%% Set-up for hpc3, 200x150x70, 200 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 10; %%% no. of x-gridpoints per tile
  sNy = 16 ; %%% no. of y-gridpoints per tile
  nPx = 20; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
   case('2km5Nodes600km')
% %  %%% Set-up for hpc3, 400x300x70, 200 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 20; %%% no. of x-gridpoints per tile
  sNy = 30 ; %%% no. of y-gridpoints per tile
  nPx = 20; %%% no. of processors in x-direction
  nPy = 10; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
     case('16km1Nodes600km')
% %  %%% Set-up for hpc3, 50x40x70, 40 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 5; %%% no. of x-gridpoints per tile
  sNy = 10 ; %%% no. of y-gridpoints per tile
  nPx = 10; %%% no. of processors in x-direction
  nPy = 4; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
  
  
    case('10km1Nodes600km')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 10; %%% no. of x-gridpoints per tile
  sNy = 12 ; %%% no. of y-gridpoints per tile
  nPx = 8; %%% no. of processors in x-direction
  nPy = 5; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
   case('20km1Nodes600km')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 5; %%% no. of x-gridpoints per tile
  sNy = 6 ; %%% no. of y-gridpoints per tile
  nPx = 8; %%% no. of processors in x-direction
  nPy = 5; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
   case('25km1Nodes600km')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 4; %%% no. of x-gridpoints per tile
  sNy = 6 ; %%% no. of y-gridpoints per tile
  nPx = 8; %%% no. of processors in x-direction
  nPy = 4; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
   case('25km1Nodes600kmNr133')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 4; %%% no. of x-gridpoints per tile
  sNy = 6 ; %%% no. of y-gridpoints per tile
  nPx = 8; %%% no. of processors in x-direction
  nPy = 4; %%% no. of processors in y-direction
  Nr = 133; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
     case('25km1Nodes600kmNr35')
% %  %%% Set-up for hpc3, 200x128x70, 160 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 4; %%% no. of x-gridpoints per tile
  sNy = 6 ; %%% no. of y-gridpoints per tile
  nPx = 8; %%% no. of processors in x-direction
  nPy = 4; %%% no. of processors in y-direction
  Nr = 35; %%% no. of z-gridpoints 
  % acct = 'cla143'; %%% XSEDE account to be charged
  cluster = 'hpc3'; %%% Destination cluster
  
  
  
  
     case('2D1km1Nodes')
% %  %%% Set-up for hpc3, 1x256x70, 40 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 1; %%% no. of x-gridpoints per tile
  sNy = 16; %%% no. of y-gridpoints per tile
  nPx = 1; %%% no. of processors in x-direction
  nPy = 32; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster 
  
     case('2D1km2Nodes')
% %  %%% Set-up for hpc3, 1x256x70, 40 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 1; %%% no. of x-gridpoints per tile
  sNy = 8; %%% no. of y-gridpoints per tile
  nPx = 1; %%% no. of processors in x-direction
  nPy = 64; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster 
  
  
     case('2D2km1Nodes')
% %  %%% Set-up for hpc3, 1x256x70, 40 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 1; %%% no. of x-gridpoints per tile
  sNy = 8; %%% no. of y-gridpoints per tile
  nPx = 1; %%% no. of processors in x-direction
  nPy = 32; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster 
  
  case('2D2km2Nodes')
% %  %%% Set-up for hpc3, 1x256x70, 80 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 1; %%% no. of x-gridpoints per tile
  sNy = 4; %%% no. of y-gridpoints per tile
  nPx = 1; %%% no. of processors in x-direction
  nPy = 64; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster 
  
   case('2D4km1Nodes')
% %  %%% Set-up for hpc3, 1x256x70, 40 procs
  opt_file = 'hpc3'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  sNx = 1; %%% no. of x-gridpoints per tile
  sNy = 4; %%% no. of y-gridpoints per tile
  nPx = 1; %%% no. of processors in x-direction
  nPy = 32; %%% no. of processors in y-direction
  Nr = 70; %%% no. of z-gridpoints 
  cluster = 'hpc3'; %%% Destination cluster 
end



    

















  %%% Uploading/downloading parameters 
  switch(cluster)
            
   case 'hpc2'
        
      username = 'hweiaf';
      clustername = 'hpc2.ust.hk';
%       toolsdir = fullfile('/home/',username,'/work/MITgcm_SM/tools/');
%       clusterdir = fullfile('/home/',username,'/work/MITgcm_SM/experiments/',batch_name); 
       toolsdir = strcat('/home/',username,'/work/MITgcm_SM/tools/');
%        clusterdir = strcat('/home/',username,'/work/MITgcm_SM/experiments/',batch_name);  
       clusterdir = strcat('/home/',username,'/work/MITgcm_SM/experiments/');    
       
 case 'hpc3'
        
      username = 'hweiaf';
      clustername = 'hpc3.ust.hk';
      toolsdir = strcat('/scratch/PI/yanwang/',username,'/work/MITgcm_SM/tools/');
      clusterdir = strcat('/scratch/PI/yanwang/',username,'/work/MITgcm_SM/experiments/');     
       
       
    otherwise %%% Defaults to Ardbeg
  
      username = 'huaiyu';
      clustername = 'caolila.atmos.ucla.edu';
      toolsdir = '/data1/MITgcm_SM/tools/';
      clusterdir = fullfile('/data1/MITgcm_SM/experiments/',batch_name);      
  
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% NON-USER-SPECIFIED PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Paths to sub-directories
%   dirpath = fullfile(runsdir,batch_name);
  dirpath = fullfile(runsdir);
  exppath = fullfile(dirpath,exp_name);
  buildpath = fullfile(exppath,builddir);
  codepath =  fullfile(exppath,codedir);
  inputpath = fullfile(exppath,inputdir);
  resultspath = fullfile(exppath,resultsdir);

  %%% We have to use MPI if we're using PBS
  if (use_pbs)
    use_mpi = true;    
  end
  
  %%% If use_mpi is false then we can only have one processor
  if ((use_mpi == false) && ((nPx~=1) || (nPy~=1)))
    error('Only one processor allowed for non-MPI computation');
  end

  %%% Calculate total grid size and number of nodes
  Nx = sNx*nSx*nPx;
  Ny = sNy*nSy*nPy;
  nodes = nPx*nPy;


  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIRECTORIES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%


  %%% Open experiment top-level directory
  [dir_success,dir_msg,dir_msgid] = mkdir(exppath);
  if (dir_success == 0)
    error(strcat(['Could not open ',exp_name,' : ',num2str(dir_msgid),' : ',dir_msg]));
  end

  %%% Open sub-directories
  subdirnames = {builddir,codedir,inputdir,resultsdir};
  for n=1:1:length(subdirnames)     
    [dir_success,dir_msg,dir_msgid] = mkdir(exppath,subdirnames{n});
    if (dir_success == 0)
      error(strcat(['Could not open ',exppath,subdirnames{n},' : ',num2str(dir_msgid),' : ',dir_msg]));
    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%


  %%% Generate 'data' and 'data.rbcs'
 eval( strcat('Nt =',setParamsNAME, '(inputpath,codepath,listterm,Nx,Ny,Nr);'));
 
 
 copyfile(strcat('./',setParamsNAME,'.m'),fullfile(inputpath));
  copyfile(strcat('./newexp_coarse.m'),fullfile(inputpath));

  %%% Generate 'eedata'
  create_eedata(inputpath,listterm,nTx,nTy);

  %%% Copy other files across
  codelist = dir('./DEFAULTS/input/');
  for n=1:1:length(codelist)
    %%% Ignore hidden files
    if (codelist(n).name(1) == '.')
      continue;
    end    
    copyfile(fullfile('./DEFAULTS/input/',codelist(n).name),fullfile(inputpath,codelist(n).name));
  end      

  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% CODE FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%


  %%% Generate SIZE.h and just copy other code files
  createSIZEh(codepath,sNx,sNy,nSx,nSy,nPx,nPy,OLx,OLy,Nr);
  codelist = dir('./DEFAULTS/code/');
  for n=1:1:length(codelist)
    if (codelist(n).name(1) == '.')
      continue;
    end
    copyfile(fullfile('./DEFAULTS/code/',codelist(n).name),fullfile(codepath,codelist(n).name));
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ESTIMATE WALL TIME %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Computation time (in hours) per gridpoint (in space and time) 
  %%% assigned to each processor.
  %%% Estimated for a single Fram core.
  alpha = 0.63e-9;
  
  %%% Estimated total computation time in hours (assumes one tile per
  %%% processor). 
  
  %%% This tends to overestimate the computation time when OLx is
  %%% comparable to sNx or OLy is comparable to sNy.
  % comptime = alpha*(sNx+2*OLx)*(sNy+2*OLy)*Nr*Nt  
  
  %%% This seems to provide a decent estimate when OLx is
  %%% comparable to sNx or OLy is comparable to sNy; 'ghost' gridpoints
  %%% require perhaps half as much processing as 'real' gridpoints.
  comptime = alpha*(sNx+OLx)*(sNy+OLy)*Nr*Nt  
    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR BUILDING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Build commands - depend on whether MPI is used
  if (use_mpi)
    mpistr = '-mpi ';
  else 
    mpistr = '';
  end
  

  

%   buildcommands = strcat([...    
%     'rm my_opt_file',lf,...
%     'ln -s ',strcat(toolsdir,'build_options/',opt_file),' my_opt_file ',lf, ...
%     strcat(toolsdir,'genmake2'),' ',mpistr,'-mods=../code -of=my_opt_file ',lf, ...
%     'make depend ',lf, ...
%     'make --always-make -j 2',lf,]);


  buildcommands = strcat([...    
    'rm my_opt_file',lf,...
    'ln -s ',strcat('../../../tools/','build_options/',opt_file),' my_opt_file ',lf, ...
    strcat('../../../tools/','genmake2'),' ',mpistr,'-mods=../code -of=my_opt_file ',lf, ...
    'make depend ',lf, ...
    'make --always-make -j 2',lf,]);


  %%% Create the 'build' shell file
  fid = fopen(fullfile(buildpath,'build.sh'),'w');
  if (fid == -1)
    error('Could not open build script for writing');
  end
  fprintf(fid,buildcommands);
  fclose(fid);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR RUNNING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Commands to link input and build folders to results folder
  runcommands = [...
    'ln -s ../',inputdir,'/* ./ ',lf, ...
    'ln -s ../',builddir,'/mitgcmuv ',lf];

  %%% Execution command depends on whether MPI is used
  if (use_pbs)
    switch (cluster)
       case 'hpc2'
        createPBSfile_hpc2(resultspath,exp_name,nodes,2*comptime,strcat(clusterdir,exp_name,'/',resultsdir));       
        runcommands = [runcommands,'sbatch run_mitgcm',lf];
        case 'hpc3'
        createPBSfile_hpc3(resultspath,exp_name,nodes,2*comptime,strcat(clusterdir,exp_name,'/',resultsdir));       
        runcommands = [runcommands,'sbatch run_mitgcm',lf];
        otherwise %%% Defaults to Ardbeg
        createPBSfile(resultspath,exp_name,nodes);
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
    end      
    
    
  else
    if (use_mpi)
      runcommands = [runcommands,'mpirun -np ',num2str(nodes), ...
                      ' ./mitgcmuv > output.txt',lf];
    else
      runcommands = [runcommands,'./mitgcmuv > output.txt',lf];
    end
  end

  %%% Create the 'run' shell file
  fid = fopen(fullfile(resultspath,'run.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,runcommands);
  fclose(fid);
  
  %%% Copy other files across
  resultslist = dir('./DEFAULTS/results/');
  for n=1:1:length(resultslist)
    %%% Ignore hidden files and run script template
    if ((resultslist(n).name(1) == '.') || strcmp(resultslist(n).name,'run_mitgcm'))
      continue;
    end    
    copyfile(fullfile('./DEFAULTS/results/',resultslist(n).name),fullfile(resultspath,resultslist(n).name));
  end    


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR COMPILING AND RUNNING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% A file to build and run MITgcm, for user convenience if they're
  %%% confident that everything is set up correctly

  %%% Commands to link input and build folders to results folder
  commands = [...
    'cd ./',builddir,'/',lf, ...
    'sh build.sh',lf, ...
    'cd ../',resultsdir,'/ ',lf, ...
    'sh run.sh',lf ];

  %%% Create the 'run' shell file
  fid = fopen(fullfile(exppath,'build_and_run.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,commands);
  fclose(fid);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILES FOR UPLOADING AND DOWNLOADING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
%   %%% Upload command
%   upcommand = [...
%     'scp -r ',...
%     '../',exp_name,'/ ', ...
%     username,'@',clustername,':',clusterdir];
%   fid = fopen(fullfile(exppath,'upload_to_cluster.sh'),'w');
%   if (fid == -1)
%     error('Could not open run script for writing');
%   end
%   fprintf(fid,upcommand);
%   fclose(fid);
%   
%   %%% Download command
%   downcommand = [...
%     'rsync -av --update ', ...
%     username,'@',clustername,':', ...
%     fullfile(clusterdir,exp_name,resultsdir),'/*ta ', ...
%     './results/ \n'];  
%   fid = fopen(fullfile(exppath,'download_from_cluster.sh'),'w');
%   if (fid == -1)
%     error('Could not open run script for writing');
%   end
%   fprintf(fid,downcommand);  
%   fclose(fid);
  

