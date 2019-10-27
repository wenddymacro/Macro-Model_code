% form P permutation from lgy_ to var_list_
P=zeros(size(lgy_,1),size(var_list_,1));
for i = 1:size(lgy_,1)
  ii = strmatch(deblank(lgy_(i,:)),var_list_,'exact');
  if ~isempty(ii)
    P(i,ii(1)) = 1;
  end
end

W = zeros(size(lgy_,1),size(lgy_,1));
W(1:size(optim_weights_,1),1:size(optim_weights_,1)) = optim_weights_;

obj = diag(W*P*oo_.var*P');
inz = find(obj > 0);
disp(' ');
disp('Objective contributions:');
disp([lgy_(inz,:) num2str(obj(inz))]);
disp(['Objective total: ' num2str(sum(obj))]);
