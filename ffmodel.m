function ret = ffmodel(NT,powr,lambda,prob,pl)

%{
'ffmodel.m' is designed to model a brush fire and output time of fire, area burned, and perimeter of burned area
Copyright (C) 2016  Samuel Cieszynski
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  
Usage: brushmodel(NT,powr,lambda,prob,pl)
  
  NT: num of time steps. If fire dies before NT, halt
  powr: grid size is 2^(powr)x2^(powr)
  lambda: expected value of a Poisson dist
  prob: prob that neighboring cells start on fire
  pl: if 1, plot matrix. If 0, don't plot
%}

%{
Outline of Program
  >Create an NxN grid of cells (Used an N^2x1 vector)
  >Assign 'fuel' to cells using a Poisson distribution with parameter lambda
  >Set random cell with fuel on fire
  >For each time step, apply the followting update rules to each cell:
    --If cell is on fire, reduce amount of 'fuel' by 1
    --If cell is on fire and has no fuel, cell is no longer on fire
    --If cell is on fire, set any 4 nearest neighbors with fuel on fire
  >Repeat previous step until no more cells are on fire
  >Record & output: time, area, perimeter, and RG calculations of A & P
  >Plot code--displays cells in a matrix
%}

N = 2^powr; % size of NxN matrix
fctr = 1; % for Plot code

state = zeros(N,N,2,2); % Matrix used to plot

% Vectors are used in update rules.
%Note: state==reshape(vector,N,N).
vector11 = zeros(N^2,1);
vector12 = zeros(N^2,1);
newvector11 = zeros(N^2,1);
newvector12 = zeros(N^2,1);

% Randomly give fuel to cells
vector12 = poissrnd(lambda,N^2,1);

vector12init = vector12;

% Start a random cell with fuel on fire
if(length(find(vector12~=0))>0)
  idx = find(vector12~=0);
  vector11(idx(randi(length(idx)))) = 1;
endif

%plot init state
if(pl)
  state(:,:,1,1) = reshape(vector11,N,N);
  state(:,:,1,2) = reshape(vector12,N,N);
  figure(1);
  hold off;
  for i =1:N
    for j = 1:N
      if (state(j,i,1,1)==0 && state(j,i,1,2)==1)
        plot(j,i,'go','MarkerSize',(18-2*powr),'MarkerFaceColor','g');hold on;title(['time=',num2str(0)]);
      endif
      if (state(j,i,1,1)==0 && state(j,i,1,2)>1)
        plot(j,i,'bo','MarkerSize',(18-2*powr),'MarkerFaceColor','b');hold on;
      endif
      if (state(j,i,1,1)~=0 && state(j,i,1,2)~=0)
        plot(j,i,'r^','MarkerSize',(15-2*powr),'MarkerFaceColor','r');hold on;
      endif
    endfor
  endfor
  axis([0 (N+1) 0 (N+1)]);drawnow
  fctr = fctr + 1;
  saveas(figure(1),['/home/samuel/Documents/octave-output/brushfire0.png'],'png');
endif

% Make sure both vectors are the same
newvector11 = vector11;
newvector12 = vector12;

% Probability of starting nearest neighbors (4 adjacent) on fire
Puunn = prob;

t = 0; % Used to measure number of time steps
while length(find(vector11==1))~=0 && t~=NT
  t = t+1;
  
  % Make boundaries not periodic
  vector11j_1 = vector11;
  vector11j_1(N:N:end) = 0;
  vector11j1 = vector11;
  vector11j1(N+1:N:end) = 0;
  
  % Burn neighboring underbrush
  newvector11([zeros(N,1);vector11(1:end-N)]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
  newvector11([vector11(1+N:end);zeros(N,1)]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
newvector11([0;vector11j_1(1:end-1)]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
newvector11([vector11j1(2:end);0]==1 & vector11==0 & vector12>0 & Puunn>rand) = 1;
  
  % For cells on fire, remove 1 fuel. If fuel==0, cell is no longer on fire
  newvector12(vector11>0 & vector12>0) -=1;
  newvector11(vector12>0 & newvector12==0) = 0;
  
  % carry over states, for looping
vector11 = newvector11;
vector12 = newvector12;
  
  if(pl) % If 1, plot the matrix at a regular interval
  state(:,:,1,1) = reshape(vector11,N,N);
  state(:,:,1,2) = reshape(vector12,N,N);
    if(mod(t,1)==0) % How often matrix is plotted
      figure(1);
      hold off;
      for i =1:N
        for j = 1:N
         if (state(i,j,1,1)==0 && state(i,j,1,2)==1)
          plot(i,j,'go','MarkerSize',(18-2*powr),'MarkerFaceColor','g');hold on;title(['time=',num2str(t)]);
        endif
        if (state(j,i,1,1)==0 && state(j,i,1,2)>1)
          plot(j,i,'bo','MarkerSize',(18-2*powr),'MarkerFaceColor','b');hold on;
        endif
        if (state(i,j,1,1)~=0 && state(i,j,1,2)~=0)
            plot(i,j,'r^','MarkerSize',(15-2*powr),'MarkerFaceColor','r');hold on;
        endif
        endfor
      endfor
    axis([0 (N+1) 0 (N+1)]);drawnow
    fctr = fctr + 1;
    endif %if(mod)
    saveas(figure(1),['/home/samuel/Documents/octave-output/brushfire', num2str(t), '.png'],'png');
  endif  %if(pl)
endwhile %end time loop


%%%%%%%%%%%%%%%%
%%  Results Calculations  %%
%%%%%%%%%%%%%%%%


Af = zeros(powr,1); % area of fire
Pf = zeros(powr,1); % perimeter of fire

bc = zeros(N^2,1); % burned cells
statebc = zeros(N,N);
Pbc = zeros(N^2,1); % perimeter of burned cells

perim = 1;
%{
There are 3 different possible definitions of perimeter:
  1. For each burned cell:
    if touching an edge, +1 for each edge
    if touching at least 1 unburned cell that is connected to an edge, +1 for each
  2. For each burned cell:
    if touching an edge, value = 1
    if touching at least 1 unburned cell that is connected to an edge value = 1
  3. For each unburned cell:
    Ideally, no burned cells are touching an edge
    if touching a burned cell, value = 1
    Separately, if burned cell is touching an edge, value = 1
For (1), count total num of edges per block.
For (2) and (3), if any cell's value==1, block = 1

The definition of perimeter used is based on perim
%}
statePbc = zeros(N,N);

bc(vector12~=vector12init) = 1;
statebc = reshape(bc,N,N);
statePbc = statebc;

%disp(statePbc);

idx = find(bc==1);
numofones = length(idx);

for n = 1:powr
  L = 2^n;
  ratio = N/L;
  
  tempA = zeros(L,L);
  idxa = zeros(L,L);
  idxalt = zeros(L,L);
  tempP = zeros(L,L);
  tempPbc = zeros(2^(2*n),1);
  
  % Find
    c = zeros(numofones,1);
    first = zeros(numofones,1);
    altfirst = zeros(numofones,1);
    second = zeros(numofones,1);
    c = mod(idx,N);
    c(c==0) = N;
    first = ceil(c/ratio);
    second = ceil((idx/N)/ratio);
    
    y = unique([first,second],"rows");
    
    i = y(:,1);
    j = y(:,2);
    
    for k = 1:length(i)
    tempA(i(k),j(k)) = 1;
    endfor
    
    %if(n==powr)
    %disp(length(find(statebc==tempA)));
    %endif
  
  tempP = tempA;
  
  tempP(:,1) += 2;
  tempP(:,end) += 2;
  tempP(1,:) += 2;
  tempP(end,:) +=2;
  tempPbc = reshape(tempP,2^(2*n),1);
  
  num1 = 1;
  num2 = 0;
  while(num1~=num2)
    num2 = num1;
    
    tempPbcj_1 = tempPbc;
    tempPbcj_1(L:L:end) = 0;
    tempPbcj1 = tempPbc;
    tempPbcj1(L+1:L:end) = 0;
    
    tempPbc([zeros(L,1);tempPbc(1:end-L)]==2 & tempPbc==0) = 2;
    tempPbc([tempPbc(1+L:end);zeros(L,1)]==2 & tempPbc==0) = 2;
    tempPbc([0;tempPbcj_1(1:end-1)]==2 & tempPbc==0) = 2;
    tempPbc([tempPbcj1(2:end);0]==2 & tempPbc==0) = 2;
    
    num1 = length(find(tempPbc==2));
  endwhile
  
  Af(n) = length(find(tempA==1));
  
  countPbc = zeros(2^(2*n),1);
  
  tempPbcj_1 = tempPbc;
  tempPbcj_1(L:L:end) = 0;
  tempPbcj1 = tempPbc;
  tempPbcj1(L+1:L:end) = 0;
  
  countPbc(tempPbc==5) = 2;
  countPbc(tempPbc==3) = 1;
  countPbc(([zeros(L,1);tempPbc(1:end-L)]==2 | [zeros(L,1);tempPbc(1:end-L)]==4) & (tempPbc==1 | tempPbc==3 | tempPbc==5)) += 1;
  countPbc(([tempPbc(1+L:end);zeros(L,1)]==2 | [tempPbc(1+L:end);zeros(L,1)]==4) & (tempPbc==1 | tempPbc==3 | tempPbc==5)) += 1;
  countPbc(([0;tempPbcj_1(1:end-1)]==2 | [0;tempPbcj_1(1:end-1)]==4) & (tempPbc==1 | tempPbc==3 | tempPbc==5)) += 1;
  countPbc(([tempPbcj1(2:end);0]==2 | [tempPbcj1(2:end);0]==4) & (tempPbc==1 | tempPbc==3 | tempPbc==5)) += 1;
  
  Pf(n) = sum(countPbc);
  
endfor

ret.t = t;
ret.Af = Af;
ret.Pf = Pf;
end %end fcn
