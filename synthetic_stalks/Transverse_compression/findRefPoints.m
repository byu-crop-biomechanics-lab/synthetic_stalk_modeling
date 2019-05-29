function [Rx1,Ry1,Rx2,Ry2] = findRefPoints(x1,y1,x2,y2,xcurve,ycurve,plotting)
% Determine the reference points that should be used in Abaqus based on the
% linear equation of the line between the points chosen on the real cross
% section.

xline = [x1 x2];
yline = [y1 y2];

% ensure that xcurve and ycurve are row vectors!

% Separate the ellipse into two parts since we're looking for two
% intercepts
xcurve1 = xcurve(1:180);
xcurve2 = xcurve(181:360);
ycurve1 = ycurve(1:180);
ycurve2 = ycurve(181:360);

% Determine xy coordinates of first reference point
b_line = polyfit(xline,yline,1);
yline2 = polyval(b_line,xcurve1);    % Evaluate 'line' at 'xcurve'

Rx1 = interp1((yline2-ycurve1),xcurve1,0);   % x-value at intercept 1
Ry1 = polyval(b_line,Rx1);                % y-value at intercept 1

% Determine xy coordinates of second reference point
% b_line = polyfit(xline,yline,1);
yline2 = polyval(b_line,xcurve2);    % Evaluate 'line' at 'xcurve'

Rx2 = interp1((yline2-ycurve2),xcurve2,0);   % x-value at intercept 1
Ry2 = polyval(b_line,Rx2);                % y-value at intercept 1

% Check with a plot
if plotting == 1
    figure(1)
    plot(xline, yline)
    hold on
    plot(xcurve1, ycurve1)
    plot(xcurve2, ycurve2)
    axis equal
    plot(Rx1, Ry1, '+r')                                % Plot Intercept Point
    plot(Rx2, Ry2, '+r')
    hold off
end

end