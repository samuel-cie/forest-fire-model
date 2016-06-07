function numoffires = modelpoiss(NT,powr,pl)

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

scale = 1; % length of a cell, in meters
ulambda = 2; % avg mass of underbrush
tlambda = 1; % avg mass of trees
Pu0 = 1; % probability of a cell having no underbrush
Pt0 = 1; % probability of a cell having no tree

% Probabilities of n starting on fire by m (Pnm)
Puunn = 1; % unn means "nearest underbrush neighbors"
Put = 0;
Putnn = 0;
Pttnn = 1;
Ptu = 0;
Ptunn = 0;

%%%%%%%%%%%%%%%%%%
%%   Dependent Variables   %%%
%%%%%%%%%%%%%%%%%%

Nc = 2^powr; % # of cells in a row
fctr = 1; % for Plot code
N = 2; % # of attributes each "cell" contains

% state is a 4D matrix. 3rd D is the distinction
% between underbrush and trees
state = zeros(Nc,Nc,2,N);
newstate = zeros(Nc,Nc,2,N);
nf = zeros(Nc*2,1);

umax = 10*scale^2; % upper bound of underbrush (in kg)
tmax = 10*scale^2; % upper bound of trees (in kg)

%%%%%%%%%%%%%%%%%%
%%%%   Initial State  %%%%%%
%%%%%%%%%%%%%%%%%%

% Provde an initial mass to all cells
state(:,:,1,2) = poissrnd(ulambda,Nc,Nc,1,1);
state(:,:,2,2) = poissrnd(tlambda,Nc,Nc,1,1);
state(Nc/2,Nc/2,1,2) = 1;

% Check if all cells are within their specified limits
if (max(state(:,:,1,2))>umax || max(state(:,:,2,2))>tmax)
  for i = 1:Nc
    for j = 1:Nc
      if (state(i,j,1,2)>umax)
        state(i,j,1,2) = umax;
      endif
      if (state(i,j,2,2)>tmax)
        state(i,j,2,2) = tmax;
      endif
    endfor
  endfor
endif

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

t = 1;
nf(1) = 1;
while nf(t)~=0 && t~=NT
t = t+1;
%for t = 1:NT
  for i =1:Nc
    for j = 1:Nc
      for k = 1:2 % 1 = underbrush; 2 = trees
        if((state(i,j,k,1)==0) && newstate(i,j,k,1)==0) %if not on fire, cells keep current value
          newstate(i,j,k,1)=state(i,j,k,1);
          newstate(i,j,k,2)=state(i,j,k,2);
          
          elseif(state(i,j,k,1)~=0 &&newstate(i,j,k,1)~=0) % burn neighbors
            
            if(k==1)
              if(state(i,j,2,1)==0 && state(i,j,2,2)~=0 && Ptu*state(i,j,k,2)>rand)
                newstate(i,j,2,1) = 1;
              endif
              if(i+1<=Nc) % If not an edge, burn
                if(state(i+1,j,1,1)==0 && state(i,j,1,2)~=0 && Puunn*state(i,j,k,2)>rand)
                  newstate(i+1,j,1,1) = 1;
                endif
                if(state(i+1,j,2,1)==0 && state(i,j,2,2)~=0 && Ptunn*state(i,j,k,2)>rand)
                  newstate(i+1,j,2,1) =1;
                endif
              endif
              if(i-1>=1) % If not an edge, burn
                if(state(i-1,j,1,1)==0 && state(i,j,1,2)~=0 && Puunn*state(i,j,k,2)>rand)
                  newstate(i-1,j,1,1) = 1;
                endif
                if(state(i-1,j,2,1)==0 && state(i,j,2,2)~=0 && Ptunn*state(i,j,k,2)>rand)
                  newstate(i-1,j,2,1) = 1;
                endif
              endif
              if(j+1<=Nc) % If not an edge, burn
                if(state(i,j+1,1,1)==0 && state(i,j,1,2)~=0 && Puunn*state(i,j,k,2)>rand)
                  newstate(i,j+1,1,1) = 1;
                endif
                if(state(i,j+1,2,1)==0 && state(i,j,2,2)~=0 && Ptunn*state(i,j,k,2)>rand)
                  newstate(i,j+1,2,1) = 1;
                endif
              endif
              if(j-1>=1) % If not an edge, burn
                if(state(i,j-1,1,1)==0 && state(i,j,1,2)~=0 && Puunn*state(i,j,k,2)>rand)
                  newstate(i,j-1,1,1) = 1;
                endif
                if(state(i,j-1,2,1)==0 && state(i,j,2,2)~=0 && Ptunn*state(i,j,k,2)>rand)
                  newstate(i,j-1,2,1) = 1;
                endif
              endif
            endif
            
            if(k==2)
              if(state(i,j,1,1)==0 && state(i,j,1,2)~=0 && Put*state(i,j,k,2)>rand)
                newstate(i,j,1,1) = 1;
              endif
              if(i+1<=Nc) % If not an edge, burn
                if(state(i+1,j,1,1)==0 && state(i,j,1,2)~=0 && Putnn*state(i,j,k,2)>rand)
                  newstate(i+1,j,1,1) = 1;
                endif
                if(state(i+1,j,2,1)==0 && state(i,j,2,2)~=0 && Pttnn*state(i,j,k,2)>rand)
                  newstate(i+1,j,2,1) = 1;
                endif
              endif
              if(i-1>=1) % If not an edge, burn
                if(state(i-1,j,1,1)==0 && state(i,j,1,2)~=0 && Putnn*state(i,j,k,2)>rand)
                  newstate(i-1,j,1,1) = 1;
                endif
                if(state(i-1,j,2,1)==0 && state(i,j,2,2)~=0 && Pttnn*state(i,j,k,2)>rand)
                  newstate(i-1,j,2,1) = 1;
                endif
              endif
              if(j+1<=Nc) % If not an edge, burn
                if(state(i,j+1,1,1)==0 && state(i,j,1,2)~=0 && Putnn*state(i,j,k,2)>rand)
                  newstate(i,j+1,1,1) = 1;
                endif
                if(state(i,j+1,2,1)==0 && state(i,j,2,2)~=0 && Pttnn*state(i,j,k,2)>rand)
                  newstate(i,j+1,2,1) = 1;
                endif
              endif
              if(j-1>=1) % If not an edge, burn
                if(state(i,j-1,1,1)==0 && state(i,j,1,2)~=0 && Putnn*state(i,j,k,2)>rand)
                  newstate(i,j-1,1,1) = 1;
                endif
                if(state(i,j-1,2,1)==0 && state(i,j,2,2)~=0 && Pttnn*state(i,j,k,2)>rand)
                  newstate(i,j-1,2,1) = 1;
                endif
              endif
            endif
            
            newstate(i,j,k,2) = state(i,j,k,2)/2;
            if(newstate(i,j,k,2)<1)
              newstate(i,j,k,2)=0;
              newstate(i,j,k,1) = 0;
            endif
        end % elseif
      end %for k
    end %for j
  end % for i
  
  if(t==2)
    %start an area of underbrush on fire
    %rstate = randi(Nc,2,1);
    %newstate(rstate(1),rstate(2),1,1) = 1;
    newstate(Nc/2,Nc/2,1,1) = 1;
    nf(t) = nf(t)+1;
    else
      %%disp(state(:,:,1,:));
      for i = 1:Nc
        for j = 1:Nc
          for k = 1:2
            if(newstate(i,j,k,1)==1 && newstate(i,j,k,2)~=0)
              nf(t) = nf(t)+1;
            endif
          endfor
      endfor
    endfor
  endif
  
  state = newstate;
  
  if(nf(t)~=0)
    disp(["nf(", num2str(t-1), ") = ", num2str(nf(t))]);
  end
  
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

%numoffires = nf;

end %end fcn
