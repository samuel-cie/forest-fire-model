function ret = callmodel(trials,lambdamin,lambdamax,numofpoints,powr,saved)

% S Cieszynski 2016
% This function is designed to loop modelpoiss.m and plot the result

prob = 1;
pl = 0;
avg = 0;
structure = 'cb';

t = lambdamin; % start at lambdamin
while t<=lambdamax % end after going through lambdamax
	each = 0;
	for i = 1:trials % call modelpoiss.m; call a part of the structure
    if(structure=='t')
      out = modelpoiss(2^(powr+1),t,powr,pl).t;
    else
      out = max(modelpoiss(2^(powr+1),t,powr,prob,pl).cb);
    endif
		each = each+out; % add up total amount
	endfor
	avg = each/trials; % find average
	figure(1) = plot(t,avg,'b'); hold on; % plot average over lambda
	t = t+(lambdamax-lambdamin)/(numofpoints-1); %loop numofpoints times
endwhile

% add title, x- and y- labels
if(structure=='t')
title (['Avg # of time steps with a  P = ', num2str(prob), ' of starting on fire']);
ylabel ("avg # of time steps");
else
title (['Avg # of cells burned with a  P = ', num2str(prob), ' of starting on fire']);
ylabel ("avg # of cells burned");
endif
xlabel ("lambda");

% save figure
saveas(figure(1),['/home/samuel/Documents/octave-output/',structure, 'P',saved,'.png'],'png');
