function area = calc_area_triangle_heron(a, b, c)
s = (a+b+c)./2;
area = sqrt(s.*(s-a).*(s-b).*(s-c));
end