function ret = modelpoiss(NT,powr,pl)

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

ulambda = 2; % avg "mass density" of underbrush (kg/m^2)
tlambda = 0; % avg "mass density" of trees (kg/m^2)
Pu0 = 1; % probability of a cell having no underbrush
Pt0 = 0; % probability of a cell having no tree

% Probabilities of n starting on fire by m (Pnm)
Puunn = 1; % unn means "nearest underbrush neighbors"
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
newstate = zeros(Nc,Nc,2,2);
nf = zeros(NT,1); %# of cells on fire
cb = zeros(NT,1); % total # of cells started on fire
cb0 = 0; % initial # of cells that are empty

%%%%%%%%%%%%%%%%%%
%%%%   Initial State  %%%%%%
%%%%%%%%%%%%%%%%%%

% Provde an initial mass to all cells
state(:,:,1,2) = poissrnd(ulambda,Nc,Nc,1,1);
state(:,:,2,2) = poissrnd(tlambda,Nc,Nc,1,1);

% Check if all cells are within their specified limits.
% Also apply probability of cell being empty.
for i = 1:Nc
    for j = 1:Nc
    % Erase some cells
      if(rand>Pu0)
        state(i,j,1,2) = 0;
      endif
      if(rand>Pt0)
        state(i,j,2,2) = 0;
      endif
      
      % Check for empty cells
      if(state(i,j,1,1)==1 || state(i,j,1,2)==0)
        cb0 = cb0+1;
      endif
      if(state(i,j,2,1)==1 || state(i,j,2,2)==0)
        cb0 = cb0+1;
      endif
    endfor
endfor

state1 = Nc/2; % randi(Nc);
state2 = Nc/2; % randi(Nc);
state(state1,state2,1,1) = 1; % start a fire at a random location
state(state1,state2,1,2) = ulambda; %Make sure cell on fire has mass

%plot init state
if(pl)
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

newstate = state; %Make all values match
t = 1; % define time. Note: t must be one ahead of the number of time steps because nf(0) does not exist.
nf(1) = 1;
cb(1) = 1;
while nf(t)~=0 && t~=NT
  t = t+1;
  for i =1:Nc % i and j represent the coordinates of an area of a 'forest'
    for j = 1:Nc
      for k = 1:2 % 1 = underbrush; 2 = trees
        if((state(i,j,k,1)==0) && newstate(i,j,k,1)==0) %if not on fire, cells keep current value
          newstate(i,j,k,2)=state(i,j,k,2);
          
        elseif(state(i,j,k,1)==1) % burn neighbors
            
            if(k==1)
              if(state(i,j,2,1)==0 && state(i,j,2,2)>0 && Ptu>rand) % Ptu*state(i,j,k,2)>rand) Can mass play a role in transition?
                newstate(i,j,2,1) = 1;
              endif
              if(i+1<=Nc) % If not an edge, burn
                if(state(i+1,j,1,1)==0 && state(i+1,j,1,2)>0 && Puunn>rand) % Puunn*state(i,j,k,2)>rand)
                  newstate(i+1,j,1,1) = 1;
                endif
                if(state(i+1,j,2,1)==0 && state(i+1,j,2,2)>0 && Ptunn>rand) % Ptunn*state(i,j,k,2)>rand)
                  newstate(i+1,j,2,1) =1;
                endif
              endif % if i+1
              if(i-1>=1) % If not an edge, burn
                if(state(i-1,j,1,1)==0 && state(i-1,j,1,2)>0 && Puunn) % Puunn*state(i,j,k,2)>rand)
                  newstate(i-1,j,1,1) = 1;
                endif
                if(state(i-1,j,2,1)==0 && state(i-1,j,2,2)>0 && Ptunn>rand) % Ptunn*state(i,j,k,2)>rand)
                  newstate(i-1,j,2,1) = 1;
                endif
              endif % if i-1
              if(j+1<=Nc) % If not an edge, burn
                if(state(i,j+1,1,1)==0 && state(i,j+1,1,2)>0 && Puunn>rand) % Puunn*state(i,j,k,2)>rand)
                  newstate(i,j+1,1,1) = 1;
                endif
                if(state(i,j+1,2,1)==0 && state(i,j+1,2,2)>0 && Ptunn) % Ptunn*state(i,j,k,2)>rand)
                  newstate(i,j+1,2,1) = 1;
                endif
              endif % if j+1
              if(j-1>=1) % If not an edge, burn
                if(state(i,j-1,1,1)==0 && state(i,j-1,1,2)>0 && Puunn>rand) % Puunn*state(i,j,k,2)>rand)
                  newstate(i,j-1,1,1) = 1;
                endif
                if(state(i,j-1,2,1)==0 && state(i,j-1,2,2)>0 && Ptunn>rand) % Ptunn*state(i,j,k,2)>rand)
                  newstate(i,j-1,2,1) = 1;
                endif
              endif % if j-1
            endif % if k==1
            
            if(k==2)
              if(state(i,j,1,1)==0 && state(i,j,1,2)>0 && Put>rand) % Put*state(i,j,k,2)>rand)
                newstate(i,j,1,1) = 1;
              endif
              if(i+1<=Nc) % If not an edge, burn
                if(state(i+1,j,1,1)==0 && state(i+1,j,1,2)>0 && Putnn>rand) % Putnn*state(i,j,k,2)>rand)
                  newstate(i+1,j,1,1) = 1;
                endif
                if(state(i+1,j,2,1)==0 && state(i+1,j,2,2)>0 && Pttnn>rand) % Pttnn*state(i,j,k,2)>rand)
                  newstate(i+1,j,2,1) = 1;
                endif
              endif % if i+1
              if(i-1>=1) % If not an edge, burn
                if(state(i-1,j,1,1)==0 && state(i-1,j,1,2)>0 && Putnn>rand) % Putnn*state(i,j,k,2)>rand)
                  newstate(i-1,j,1,1) = 1;
                endif
                if(state(i-1,j,2,1)==0 && state(i-1,j,2,2)>0 && Pttnn>rand) % Pttnn*state(i,j,k,2)>rand)
                  newstate(i-1,j,2,1) = 1;
                endif
              endif % if i-1
              if(j+1<=Nc) % If not an edge, burn
                if(state(i,j+1,1,1)==0 && state(i,j+1,1,2)>0 && Putnn>rand) % Putnn*state(i,j,k,2)>rand)
                  newstate(i,j+1,1,1) = 1;
                endif
                if(state(i,j+1,2,1)==0 && state(i,j+1,2,2)>0 && Pttnn>rand) % Pttnn*state(i,j,k,2)>rand)
                  newstate(i,j+1,2,1) = 1;
                endif
              endif % if j+1
              if(j-1>=1) % If not an edge, burn
                if(state(i,j-1,1,1)==0 && state(i,j-1,1,2)>0 && Putnn>rand) % Putnn*state(i,j,k,2)>rand)
                  newstate(i,j-1,1,1) = 1;
                endif
                if(state(i,j-1,2,1)==0 && state(i,j-1,2,2)>0 && Pttnn>rand) % Pttnn*state(i,j,k,2)>rand)
                  newstate(i,j-1,2,1) = 1;
                endif
              endif % if j-1
            endif % if k==2
            
            newstate(i,j,k,2) = state(i,j,k,2)-1;
            if(newstate(i,j,k,2)<1)
              newstate(i,j,k,1) = 0;
            endif
        end % elseif
      end %for k
    end %for j
  end % for i
  
  for i = 1:Nc % Loop through all cells
    for j = 1:Nc
      for k = 1:2
        if(newstate(i,j,k,1)==1 && newstate(i,j,k,2)~=0)
          nf(t) = nf(t)+1; % count # of cells on fire
        endif
        if(newstate(i,j,k,1)==1 || newstate(i,j,k,2)==0)
          cb(t) = cb(t)+1; % count total # of cells burned
        endif
      endfor
    endfor
  endfor
  
  % carry over states, for looping
  state = newstate;
  
  cb(t) = cb(t)-cb0;
  
  if(pl)
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
ret.numoffires = nf;
ret.burnedcells = cb;

end %end fcn
