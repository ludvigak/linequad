function zr = rotate_and_scale(za, zb, z)
% zr = rotate_and_scale(za, zb, z)
%
% zr = M(z)
% where
% M(za) = -1
% M(zb) = 1

c_map = [ (za+zb)/(za-zb) 2/(zb-za) ];
zr = c_map(1) + z*c_map(2);
