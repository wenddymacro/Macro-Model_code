d = dbload('data.csv');

varlist={'RS_US','LCPI_US','LGDP_US','UNR_US', 'BLT_US'};
rng = qq(1994,1):qq(2007,4);

fid = fopen('data.m','w');
for i = 1:length(varlist)
  fprintf(fid, '%s = [\n', varlist{i});
  fprintf(fid, '%16.12g\n', d.(varlist{i})(rng));
  fprintf(fid, '];\n\n');
end

fclose(fid);