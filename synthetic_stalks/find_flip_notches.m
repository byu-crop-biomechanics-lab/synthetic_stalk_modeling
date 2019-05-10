clear;
close;

load Section_slices_bottom_990.mat

R_ext = squeeze(ext_rhoDCSR);
R_int = squeeze(int_rhoDCSR);
T = tDCSR';

N = size(R_ext,2);
flip = zeros(N,1);

for i = 1:N
    % Plot each cross section to see if it needs to be flipped 180 degrees
    polarplot(T,R_ext(:,i));
    i
    s = input('Enter 1 if cross section needs to flip: ');
    if isempty(s)
        s = 0;
    end
    flip(i) = s;
%     pause(); 
end