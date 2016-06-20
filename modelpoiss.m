function ret = modelpoiss(NT,ulambda,powr,prob,pl)

%S Cieszynski 2016
%Based on bakfire.m by T Skorczewski 2016
%Plot code was created by T Skorczewski 2016

% This program is a CA model of a forest fire.
% There are 2 layers in the forest: trees, and undergrowth.
% This starts a fire on the ground and lets the fire run until it dies out.

%%%%%%%%%%%%%%%%%%
%%   Independent Variables   %%
%%%%%%%%%%%%%%%%%%

% NT - number of time steps
% If NT~=0, run until NT-1. If Nt==0, run until fire burns out.
% powr - log_2(length of the matrix)
% pl - used to turn plotting on or off
% output: numoffires - # of fires located in grid

%ulambda = 2; % avg "mass density" of underbrush (kg/m^2)
tlambda = 0; % avg "mass density" of trees (kg/m^2)

% Probabilities of n starting on fire by m (Pnm)
Puunn = prob; % unn means "nearest underbrush neighbors"
Put = 0;
Putnn = 0;
Pttnn = 0;
Ptu = 0;
Ptunn = 0;

%%%%%%%%%%%%%%%%%%
%%   Dependent Variables   %%%
%%%%%%%%%%%%%%%%%%

NT = NT+1; % this is done in response to nf(0) being undefined.
Nc = 2^powr; % # of cells in a row
fctr = 1; % for Plot code

% state & newstate are both a 4D matrix.
% 1D &2D represent surface areas where a fire can occur.
% 3D represents whether a cell is underbrush (1) or a tree (2).
% 4D stores the values for whether on fire (1) and mass of cell (2).
state = zeros(Nc,Nc,2,2);

vector11 = zeros(Nc^2,1);
vector12 = zeros(Nc^2,1);
vector21 = zeros(Nc^2,1);
vector22 = zeros(Nc^2,1);
newvector11 = zeros(Nc^2,1);
newvector12 = zeros(Nc^2,1);
newvector21 = zeros(Nc^2,1);
newvector22 = zeros(Nc^2,1);

%disp(vector12);

nf = zeros(NT,1); %# of cells on fire
cb = zeros(NT,1); % total # of cells started on fire


%%%%%%%%%%%%%%%%%%
%%%%   Initial State  %%%%%%
%%%%%%%%%%%%%%%%%%

% Provde an initial mass to all cells
vector12 = poissrnd(ulambda,Nc^2,1);
vector22 = poissrnd(tlambda,Nc^2,1);

% Start a random cell on fire
if(length(find(vector12~=0))>0)
  idx = find(vector12~=0);
  vector11(idx(randi(length(idx)))) = 1;
endif

%plot init state
if(pl)
  state(:,:,1,1) = reshape(vector11,Nc,Nc);
  state(:,:,1,2) = reshape(vector12,Nc,Nc);
  figure(1);
  hold off;
  for i =1:Nc
    for j = 1:Nc
      if (state(i,j,1,1)==0 && state(i,j,1,2) ~=0)
        plot(i,j,"ms","MarkerSize",(18-2*powr),"MarkerFaceColor","m");hold on;
      end % if state(1)==0
      if(state(i,j,2,1)==0 && state(i,j,2,2) ~=0)
        plot(i,j,"go","MarkerSize",(15-2*powr),"MarkerFaceColor","g");hold on;
      end % if state(2)==0
      if (state(i,j,1,1)~=0 && state(i,j,1,2) ~=0)
        plot(i,j,"r^","MarkerSize",(15-2*powr),"MarkerFaceColor","r");hold on;
      end % if state(1)~=1
      if (state(i,j,2,1)~=0 && state(i,j,2,2) ~=0)
        plot(i,j,"k^","MarkerSize",(12-2*powr),"MarkerFaceColor","k");hold on;
      end % if state(2)~=1
    end
  end
  axis([0 (Nc+1) 0 (Nc+1)]);drawnow
  fctr = fctr + 1;
end
%%%%%%%%%%%

% Make sure both vectors are the same
newvector11 = vector11;
newvector12 = vector12;
cb0 = length(find((newvector11==1 | newvector12==0)))-1; % initial # of cells that are empty
t = 1; % define time. Note: t must be one ahead of the number of time steps because nf(0) does not exist.
nf(1) = 1;
cb(1) = 1;
while nf(t)~=0 && t~=NT
  t = t+1;
  
  % Make boundaries not periodic
  vector11j_1 = vector11;
  vector11j_1(Nc:Nc:end) = 0;
  vector11j1 = vector11;
  vector11j1(Nc+1:Nc:end) = 0;
  
  vector21j_1 = vector21;
  vector21j_1(Nc:Nc:end) = 0;
  vector21j1 = vector21;
  vector21j1(Nc+1:Nc:end) = 0;
  
  %% Burn Underbrush (directions correspond to plot)
  % burn cell above
  newvector11([zeros(Nc,1);vector11(1:end-Nc)]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
  % burn cell below
  newvector11([vector11(1+Nc:end);zeros(Nc,1)]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
  % burn cell to the left
	newvector11([0;vector11j_1(1:end-1)]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
  % burn cell to the right
	newvector11([vector11j1(2:end);0]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
  
  %% Burn Trees (directions correspond to plot)
  % burn cell above
  newvector21([zeros(Nc,1);vector21(1:end-Nc)]==1 & vector21==0 & vector22>0 & Puunn>rand) = 1;
  % burn cell below
  newvector21([vector21(1+Nc:end);zeros(Nc,1)]==1 & vector21==0 & vector22>0 & Puunn>rand) = 1;
  % burn cell to the left
	newvector21([0;vector21j_1(1:end-1)]==1 & vector21==0 & vector22>0 & Puunn>rand) = 1;
  % burn cell to the right
	newvector21([vector21j1(2:end);0]==1 & vector21==0 & vector22>0 & Puunn>rand) = 1;
  
  newvector12(vector11>0 & vector12>0) -=1;
  newvector11(vector12>0 & newvector12==0) = 0;
  
  newvector22(vector21>0 & vector22>0) -=1;
  newvector21(vector22>0 & newvector22==0) = 0;
  
  nf(t) = length(find((newvector11==1 & newvector12~=0)));
  cb(t) = length(find((newvector11==1 | newvector12==0)))-cb0;
  
  % carry over states, for looping
	vector11 = newvector11;
	vector12 = newvector12;
  vector21 = newvector21;
  vector22 = newvector22;
  
  if(pl)
  state(:,:,1,1) = reshape(vector11,Nc,Nc);
  state(:,:,1,2) = reshape(vector12,Nc,Nc);
    if(mod(t,1)==0)    
      figure(1);
      hold off;
      for i =1:Nc
        for j = 1:Nc
         if (state(i,j,1,1)==0 && state(i,j,1,2) ~=0)
          plot(i,j,'ms','MarkerSize',(18-2*powr),'MarkerFaceColor','m');hold on;title(['time=',num2str(t-1)]);
        end % if
        if(state(i,j,2,1)==0 && state(i,j,2,2) ~=0)
            plot(i,j,'go','MarkerSize',(15-2*powr),'MarkerFaceColor','g');hold on;
        end % if
        if (state(i,j,1,1)~=0 && state(i,j,1,2) ~=0)
            plot(i,j,'r^','MarkerSize',(15-2*powr),'MarkerFaceColor','r');hold on;
        end % if
        if (state(i,j,2,1)~=0 && state(i,j,2,2) ~=0)
            plot(i,j,'k^','MarkerSize',(12-2*powr),'MarkerFaceColor','k');hold on;
        end % if
        end % for j
      end % for i
    axis([0 (Nc+1) 0 (Nc+1)]);drawnow
    fctr = fctr + 1;
    end %if(mod)
  end  %if(pl)
end %end time loop

ret.t = t-1;
% first value is nf(0), or the initial conditions.
ret.nf = nf;
ret.cb = cb;

end %end fcn
