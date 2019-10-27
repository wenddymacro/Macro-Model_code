function [] = resize_subplot3(h1, h2, h3, h4, h5, h6, h7, h8, h9, rat1, rat2)
  p1 = get(h1, 'Position');
  p2 = get(h2, 'Position');
  p3 = get(h3, 'Position');
  
  p4 = get(h4, 'Position');
  p5 = get(h5, 'Position');
  p6 = get(h6, 'Position');
  
  p7 = get(h7, 'Position');
  p8 = get(h8, 'Position');
  p9 = get(h9, 'Position');
  
  tot_height = p1(2)+p4(2)+p7(2);
  height = tot_height*4/5;
  space = (tot_height-height)/5;
  p3(2) = 2*space;
  p3(4) = (1-rat1-rat2)*height;
  p2(2) = p3(2) + p3(4) + space;
  p2(4) = rat2*height;
  p1(2) = p2(2) + p2(4) + space;
  p1(4) = rat1*height;
  
  set(h1,'Position',p1);
  set(h2,'Position',p2);
  set(h3,'Position',p3);
  return;
  
