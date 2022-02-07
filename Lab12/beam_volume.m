function V = beam_volume(x)

global L

t = x(1);
w = x(2);
V = L*t*w;

end