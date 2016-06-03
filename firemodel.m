function numoffires = firemodel(NT,powr,pl)

%S Cieszynski 2016
%Based on bakfire.m by T Skorczewski 2016
%Plot code was created by T Skorczewski 2016

% This program is a CA model of a forest fire

% NT - number of time steps
% powr - 2^length of the matrix
% pl - used to turn plotting on or off
% output: numoffires - # of fires located in grid

Nc = 2^powr; % # of cells in a row
fctr = 1; % for Plot code
N = 2; % # of attributes each "cell" contains

% state is a 4D matrix. 3rd D is the distinction
% between underbrush and trees
state = zeros(Nc,Nc,2,N);
newstate = zeros(Nc,Nc,2,N);
nf = zeros(NT,1);
scale = 10; % length of a cell, in meters

% Probabilities of n starting on fire by m (Pnm)
Puunn = 1; % unn means "nearest neighbor underbrush"
Put = 0;
Putnn = 0.1;

Pttnn = 1;
Ptu = 0.2;
Ptunn = 0;

umax = 20*scale^2; % upper bound of underbrush (in kg)
tmax = 200*scale^2; % upper bound of trees (in kg)

freq = 100; % # of time steps in between start of fire
growth = 100; % # of time steps between growth of forest
growthrate = 2;
Pgrowth = 1;

burnt = 3; % # of time steps a tree is on fire
burnu = 2; % # of time steps underbrush is on fire

P0u = 0.8;
P0t = 0.8;
for i = 1:Nc %plot intial forest
  for j = 1:Nc
      if(P0u>rand)
        state(i,j,1,2) = randi(umax);
      endif
      if(P0t>rand)
        state(i,j,2,2) = randi(tmax);
      endif
  endfor
endfor

%plot init state
if(pl)
figure(1);
hold off;
for i =1:Nc
    for j = 1:Nc
        if (state(i,j,1,1)==0 && state(i,j,1,2) ~=0)
          plot(i,j,'ms','MarkerSize',(18-2*powr),'MarkerFaceColor','m');hold on;
        end % if state(1)==0
        if(state(i,j,2,1)==0 && state(i,j,2,2) ~=0)
            plot(i,j,'go','MarkerSize',(15-2*powr),'MarkerFaceColor','g');hold on;
        end % if state(2)==0
        if (state(i,j,1,1)~=0 && state(i,j,1,2) ~=0)
            plot(i,j,'r^','MarkerSize',(15-2*powr),'MarkerFaceColor','r');hold on;
        end % if state(1)~=1
        if (state(i,j,2,1)~=0 && state(i,j,2,2) ~=0)
            plot(i,j,'k^','MarkerSize',(12-2*powr),'MarkerFaceColor','k');hold on;
        end % if state(2)~=1
    end
end
axis([0 (Nc+1) 0 (Nc+1)]);drawnow
fctr = fctr + 1;
end
%%%%%%%%%%%

for t = 1:NT
  for i =1:Nc
    for j = 1:Nc
      for k = 1:2 % 1 = underbrush; 2 = trees
        if((state(i,j,k,1)==0) && newstate(i,j,k,1)==0) %if not on fire, cells keep current value
          newstate(i,j,k,1)=state(i,j,k,1);
          newstate(i,j,k,2)=state(i,j,k,2);
          
          elseif(state(i,j,k,1)~=0 &&newstate(i,j,k,1)~=0) % burn neighbors
            
            if(k==1)
              if(state(i,j,2,1)==0 && state(i,j,2,2)~=0 && Ptu>rand)
                newstate(i,j,2,1) = burnt;
              endif
              if(i+1<=Nc) % If not an edge, burn
                if(state(i+1,j,1,1)==0 && state(i,j,1,2)~=0 && Puunn>rand)
                  newstate(i+1,j,1,1) = burnu;
                endif
                if(state(i+1,j,2,1)==0 && state(i,j,2,2)~=0 && Ptunn>rand)
                  newstate(i+1,j,2,1) = burnt;
                endif
              endif
              if(i-1>=1) % If not an edge, burn
                if(state(i-1,j,1,1)==0 && state(i,j,1,2)~=0 && Puunn>rand)
                  newstate(i-1,j,1,1) = burnu;
                endif
                if(state(i-1,j,2,1)==0 && state(i,j,2,2)~=0 && Ptunn>rand)
                  newstate(i-1,j,2,1) = burnt;
                endif
              endif
              if(j+1<=Nc) % If not an edge, burn
                if(state(i,j+1,1,1)==0 && state(i,j,1,2)~=0 && Puunn>rand)
                  newstate(i,j+1,1,1) = burnu;
                endif
                if(state(i,j+1,2,1)==0 && state(i,j,2,2)~=0 && Ptunn>rand)
                  newstate(i,j+1,2,1) = burnt;
                endif
              endif
              if(j-1>=1) % If not an edge, burn
                if(state(i,j-1,1,1)==0 && state(i,j,1,2)~=0 && Puunn>rand)
                  newstate(i,j-1,1,1) = burnu;
                endif
                if(state(i,j-1,2,1)==0 && state(i,j,2,2)~=0 && Ptunn>rand)
                  newstate(i,j-1,2,1) = burnt;
                endif
              endif
            endif
            
            if(k==2)
              if(state(i,j,1,1)==0 && state(i,j,1,2)~=0 && Put>rand)
                newstate(i,j,1,1) = burnu;
              endif
              if(i+1<=Nc) % If not an edge, burn
                if(state(i+1,j,1,1)==0 && state(i,j,1,2)~=0 && Putnn>rand)
                  newstate(i+1,j,1,1) = burnu;
                endif
                if(state(i+1,j,2,1)==0 && state(i,j,2,2)~=0 && Pttnn>rand)
                  newstate(i+1,j,2,1) = burnt;
                endif
              endif
              if(i-1>=1) % If not an edge, burn
                if(state(i-1,j,1,1)==0 && state(i,j,1,2)~=0 && Putnn>rand)
                  newstate(i-1,j,1,1) = burnu;
                endif
                if(state(i-1,j,2,1)==0 && state(i,j,2,2)~=0 && Pttnn>rand)
                  newstate(i-1,j,2,1) = burnt;
                endif
              endif
              if(j+1<=Nc) % If not an edge, burn
                if(state(i,j+1,1,1)==0 && state(i,j,1,2)~=0 && Putnn>rand)
                  newstate(i,j+1,1,1) = burnu;
                endif
                if(state(i,j+1,2,1)==0 && state(i,j,2,2)~=0 && Pttnn>rand)
                  newstate(i,j+1,2,1) = burnt;
                endif
              endif
              if(j-1>=1) % If not an edge, burn
                if(state(i,j-1,1,1)==0 && state(i,j,1,2)~=0 && Putnn>rand)
                  newstate(i,j-1,1,1) = burnu;
                endif
                if(state(i,j-1,2,1)==0 && state(i,j,2,2)~=0 && Pttnn>rand)
                  newstate(i,j-1,2,1) = burnt;
                endif
              endif
            endif
            
            newstate(i,j,k,1)=state(i,j,k,1)-1;
            if(newstate(i,j,k,1)==0)
              newstate(i,j,k,2)=0;
            endif
        end % elseif
      end %for k
    end %for j
  end % for i
  
  if(mod(t,growth)==0) % All cold cells grow
    for i = 1:Nc
      for j = 1:Nc
        for k = 1:2
          newstate(i,j,k,2) = state(i,j,k,2)+growthrate*scale^2;
          if(k==1 && newstate(i,j,k,2)>umax)
            newstate(i,j,k,2) = umax;
            elseif(newstate(i,j,k,2)>tmax)
            newstate(i,j,k,2) = tmax;
          endif
        endfor %for k
      end % for j
    end % for i
  endif % if(mod)
  
  if(mod(t,freq)==1) % Start a tree on fire at random location
    rstate = randi(Nc,2,1);
    newstate(rstate(1),rstate(2),1,1) = 1;
    newstate(rstate(1),rstate(2),2,1) = 1;
  endif
  
  state = newstate;
  nf(t) = length(find(state==1));
  
  if(pl)
    if(mod(t,1)==0)    
      figure(1);
      hold off;
      for i =1:Nc
        for j = 1:Nc
         if (state(i,j,1,1)==0 && state(i,j,1,2) ~=0)
          plot(i,j,'ms','MarkerSize',(18-2*powr),'MarkerFaceColor','m');hold on;title(['time=',num2str(t)]);
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
    %mov(fctr) = getframe(gcf);
    fctr = fctr + 1;
    end %if(mod)
  end  %if(pl)
end %end time loop

numoffires = nf;

end %end fcn