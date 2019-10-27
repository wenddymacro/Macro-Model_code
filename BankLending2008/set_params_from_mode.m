% set parameters from mode given in the estim_oo_file

estim_oo_file = '../' estim_oo_file '.mat';
a = load(estim_oo_file);
estp = a.oo_.posterior_mode.parameters;
ests = a.oo_.posterior_mode.shocks_std;

pnames = fieldnames(estp);
for i = 1:length(pnames)
  evalin('base', [pnames{i} '=estp.(' pnames{i} ')']);
end

