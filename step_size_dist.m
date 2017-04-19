% fid = fopen('list_corr.txt', 'r');
close all;
clear all;
tdfread('list_corr.txt');

t_scale = 1/60; %60 fps

n_files = length(workstate_000800x2Edat);
step_size = cell(n_files,1);
msd = cell(n_files,1);
del_t = [1:1:24];


% count = 1;
max_time = [];
for i = 1 : n_files
  file = workstate_000800x2Edat(i,:);
  M = load(file);
  pos = [M(:,2),M(:,3)];
  % pause
  msd{i} = zeros(size(pos,1)-1,1);
  max_time = [max_time size(pos,1)-1];
  for j = 2 : size(pos,1)
    msd{i}(j-1) = sqrt((pos(j,1)-pos(1,1))^2 + (pos(j,2)-pos(1,2))^2);
  end
  step_size{i} = zeros(size(pos,1)-1,length(del_t));
  for k = 1 : length(del_t)
    for j = 1+del_t(k) : size(pos,1)
      step_size{i}(j-del_t(k),del_t(k)) = sqrt((pos(j,1)-pos(j-del_t(k),1))^2 + (pos(j,2)-pos(j-del_t(k),2))^2);
    end
  end
    % count = count + 1;
end


% cum_step_size = ;
% cum_msd = [];
binfellows = [0.05:0.05:8];
cols = jet(length(del_t));
figure;
msdvals = zeros(size(del_t));
for k = 1 : length(del_t)
  cum_step_size = [];
  for i = 1 : n_files
    cum_step_size = [cum_step_size; step_size{i}(1:size(step_size{i},1)-del_t(k), del_t(k))];
  end
  msdvals(k) = mean(cum_step_size.*cum_step_size);
  cum_step_size = cum_step_size/sqrt(del_t(k)^(.93));
  % figure(1)
  [f x] = hist(cum_step_size,30);
  plot(x, f/trapz(x,f),'-','Color',cols(k,:))
  % figrue(2)
  % [f,x] = ecdf(cum_step_size);
  % ones(length(f),1) - f;
  % plot(x,f,'-')
  set(gca, 'Yscale', 'log')
  % onemincdf = zeros(size(binfellows));
  % for binno = 1:length(binfellows)
  %   onemincdf(binno) = 1-(length(find(cum_step_size<=binfellows(binno)))/length(cum_step_size));
  % end
  %
  % semilogy(binfellows,onemincdf,'o','Color',cols(k,:));
  hold on
  % pause
end
dval = .48;
% semilogy(binfellows,binfellows.*exp(-(binfellows.^2)/(4*dval))/(2*dval),'-k','LineWidth',2)
hold off
axis([0 8 0.0001 1.2])
colorbar
% legend('1)

figure
loglog(del_t,msdvals,'rs')


% cum_msd = [cum_msd; msd{i}];
% plot([1:1:max(max_time)],msd{i},'*-')
% hold on
% pause
