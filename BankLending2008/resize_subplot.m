function [] = resize_subplot(h1, h2, upper_rat)
  p1 = get(h1, 'Position');
  p2 = get(h2,'Position');
  steal = upper_rat*(p1(4)+p2(4))-p1(4);
  set(h1,'Position',[p1(1) p1(2)-steal p1(3) p1(4)+steal]);
  set(h2,'Position',[p2(1:3) p2(4)-steal]);
  return;
  
