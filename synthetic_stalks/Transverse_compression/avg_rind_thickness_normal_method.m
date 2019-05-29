function [avgrindthickness, int_X, int_Y] = avg_rind_thickness_normal_method(I, ext_X, ext_Y, plotting)

% This function calculates the average rind thickness of an individual slice (I). 
% The thickness is calculated using a series of points that march inward from the exterior boundary.
% Inward marching continues until the median value of the interior curve matches the median valule of 
% the exterior curve. Along the way, wmoothing and point elimination/replacement are used to 
% insure that the interior curves are a good approximation of the original exterior shape. 
%
% INPUTS:           I - the image slice to be analyzed
%                   ext_X and ext_Y - exterior XY contour of I
%                   plotting - plotting on(1) or off(0)
%
% OUTPUTS:      avgrindthickness - the estimated average rind thickness
%               int_X and int_Y  - the interior and exterior boundary curves
%
% 2017.09.19 - No function changes. However, I added additional comments to make 
% the method readable. Major sections are now readily apparent.


% ================= SETTINGS ====================================================
% these values should be consistent.
inc = -1;           % increment to step inward
dist_in = 150;       % distance (in pixels) to step inward from exterior contour.
skip = 0;           % number of exterior points to skip before starting. 
mingap = 2;         % minimum gap between adjacent points (pixels)
midgap = sqrt(8);   % medium gap between adjacent points (pixels)
maxgap = 4;         % maximum gap between adjacent points (pixels)
span = 9;           % smoothing span
minnum = 10;        % minimum number of contour points required for analysis (exit loop if fewer are available).
% ================= SETTINGS ====================================================

% INITIALIZE VARIABLES
avgrind = 0;        % average rind thickness variable.
crossing = 0;       % variable indicating whether or not the median rind thickness has crossed the medrindref value.
medrind = 0;        % variable that contains the median rind intensity at a given inset from the exterior


% ================= PRELIMINARIES ===============================================
% SKIPPING SOME DATA POINTS (OPTIONAL)
ext_X2 = ext_X(1:skip+1:end);        % keep only every third to avoid duplication of sample points.
ext_Y2 = ext_Y(1:skip+1:end);

% INITIAL SMOOTHING
S = 3; % smoothing scale factor (more smoothing at the beginning helps avoid curve overlap).
N = length(ext_X2);
if N > minnum  % check that length of ext_X2 is at least minnum before smoothing.
    hspan = (span-1)/2;             % smoothing half span

    % SMOOTHING THE ORIGINAL DATA
    ext_X2 = smooth(ext_X2(1:1:end),span*3,'sgolay',2);  % smooth external curves so that the normal directions are not affected by noise.
    ext_Y2 = smooth(ext_Y2(1:1:end),span*3,'sgolay',2);
    lapx = smooth([ext_X2(end-(span-1):end);ext_X2(1:span)],span,'sgolay',2);   % overlapping section consisting of end and begining
    lapy = smooth([ext_Y2(end-(span-1):end);ext_Y2(1:span)],span,'sgolay',2);
    ext_X2([N-hspan+1:N,1:hspan]) = lapx(hspan+2:span+hspan);   % extracting sections from the lap data to be used in ext_X2,Y2 to smooth the begining and end.
    ext_Y2([N-hspan+1:N,1:hspan]) = lapy(hspan+2:span+hspan);   
end
% ================= PRELIMINARIES ===============================================


for i = 1:dist_in  % loops inward until 

    % ===================  STEP 1 ==============================================================================
    % FILL IN ANY GAPS WHERE SAMPLING IS TOO SPARSE 
    r = sqrt(diff(ext_X2).^2 + diff(ext_Y2).^2);                                % distances between adjacent contour points
    Rfirstlast = sqrt((ext_X2(1)-ext_X2(end))^2 + (ext_Y2(1)-ext_Y2(end))^2);   % distance between first and alast points.
    r = [r; Rfirstlast];                                                        % r appended for complete list of distances between adjacent points.
    INDEX = r > maxgap;                                                         % logical index in r, ext_X2, ext_Y2 system. (this method faster than find()
    list = 1:(length(r)+1);                                                     % list of possible indices
    INDEX = list(INDEX);                                                        % actual indices of distances greater than 2

    % This if statement is only entered if gaps are identified.
    if ~isempty(INDEX)                  % if INDEX is empty, none of the following is necessary

        difs = diff(INDEX);             % differences between neighboring indices
        gaps = find(difs>2);            % indices of INDEX indicating breaks between seqential groups of indices
        starts = [1, gaps+1];           % the starting indices of INDEX for each group
        finis = [gaps, length(INDEX)];  % the finishing indices of INDEX for each group.

        for k = 1:length(starts);       % loop through the number of groups 
            if  starts(k) == finis(k) && INDEX(starts(k))~=length(ext_X2) % skip cases where there is just a single isolated gap that is not the end/begining interface
            continue, end 

            rindex = INDEX(starts(k)):INDEX(finis(k)); % indices of ext_X2, r, and ext_Y2 in each group. Note the nested indices.   
            xindex = [rindex, rindex(end)+1];          % xindex is the indices of ext_X2,Y2, it is one longer than rindex because rindex is based on differences between adjacent values
            if xindex(end) > length(ext_X2)            % of xindex are greater than the size of xindex, these refer to the first index. 
                xindex(end) = 1;
            end

            t = [0; cumsum(r(rindex))];                 % cumulative sum of distances (pathwise coordinate system)
            total = t(end);                             % the total length between start and finish
            n = round(total/midgap);                    % number of new points 
            if n == 1 | n == 2, n = n + 1; end          % if n is 1 or 2, increase n by 1.

            x = ext_X2(xindex);                         % data points
            y = ext_Y2(xindex);                         % ditto
            sample = linspace(t(1),t(end),n);           % new sample points
            xnew = interp1(t,x,sample);                 % new data points
            ynew = interp1(t,y,sample);                 % ditto

            if xindex(end) == 1                         % the case of beginning/end interface
                ext_X2 = [ext_X2(1:xindex(1)); xnew(2:end-1)'];             % insert new x values
                ext_Y2 = [ext_Y2(1:xindex(1)); ynew(2:end-1)'];             % ditto
                r = [r(1:xindex(1)); 0*ynew(2:end-1)'; r(xindex(end):end)]; % zeros inserted to adjust r to align with ext_X2 and ext_Y2
            else
                ext_X2 = [ext_X2(1:xindex(1)); xnew(2:end-1)'; ext_X2(xindex(end):end)];    % insert new x values
                ext_Y2 = [ext_Y2(1:xindex(1)); ynew(2:end-1)'; ext_Y2(xindex(end):end)];    % ditto
                r = [r(1:xindex(1)); 0*ynew(2:end-1)'; r(xindex(end):end)];                 % zeros inserted to adjust r to align with ext_X2 and ext_Y2
            end

            INDEX = INDEX + n - length(xindex);         % adjust index values to account for increases in the length of ext_X2,Y2.

        end
    end
    % ===================  END STEP 1 ==============================================================================



    % ===================  STEP 2 ==============================================================================
    % ==== REMOVE ANY THAT ARE TOO CLOSE TOGETHER  =======================================================
    r = sqrt(diff(ext_X2).^2 + diff(ext_Y2).^2);    % vector of distances between adjacent points.
    N = length(ext_X2);     % length
    kill = zeros(1,N);      % index of which points to "kill" (remove).
    k = 1;                  % start at k = 1 
    while k < N             % look through all k values      
        if r(k) > mingap    % if the distances is greater than the min, just keep going 
            k = k + 1;      % increment k
        else
            kill(k+1) = 1;  % else (the gap is smaller/= than mingap)
            j = 2;          % initialize j
            if k + j > N    % don't proceed past k + j = N     
                break;      
            else
                while sqrt((ext_X2(k)-ext_X2(k+j))^2 + (ext_Y2(k)-ext_Y2(k+j))^2) < mingap  % if/while the distance is less than mingap
                    kill(k+j) = 1;                  % record the index of offending points
                    j = j + 1;                      % increment j
                    if k + j > N, break; end        % don't proceed beyond k + j = N
                end
            end
            k = k + j;  % update k to current index
        end
    end

    keep = logical(1 - kill);   % invert kill indicess to create "keep" indices
    Rfirstlast = sqrt((ext_X2(1)-ext_X2(N))^2 + (ext_Y2(1)-ext_Y2(N))^2);   % calculate distance between first and last points
    if Rfirstlast < mingap      % don't keep if it's below mingap
        keep(N) = 0;
    end

    ext_X2 = ext_X2(keep);      % keep all values below the threshold      
    ext_Y2 = ext_Y2(keep);      % ditto   
    % ===================  END STEP 2 ==============================================================================


    % ===================  STEP 3 ==============================================================================
    % SMOOTH
    N = length(ext_X2);
    if N > minnum  % check that length of ext_X2 is at least minnum before smoothing.
        hspan = (span-1)/2;             % smoothing half span
        ext_X2 = smooth(ext_X2(1:1:end),span,'sgolay',2);  % smooth external curves so that the normal directions are not affected by noise.
        ext_Y2 = smooth(ext_Y2(1:1:end),span,'sgolay',2);
        lapx = smooth([ext_X2(end-(span-1):end);ext_X2(1:span)],span,'sgolay',2);   % overlapping section consisting of end and begining
        lapy = smooth([ext_Y2(end-(span-1):end);ext_Y2(1:span)],span,'sgolay',2);
        ext_X2([N-hspan+1:N,1:hspan]) = lapx(hspan+2:span+hspan);   % extracting sections from the lap data to be used in ext_X2,Y2 to smooth the begining and end.
        ext_Y2([N-hspan+1:N,1:hspan]) = lapy(hspan+2:span+hspan);   
    end
    % ===================  END STEP 3 ==============================================================================



    % ===================  STEP 4 ==============================================================================
    % NORMAL ANGLES
    normal_angle = zeros(N,1);              % initialize variable
    normal_angle(2:N-1) = atan2(-(ext_X2(3:N) - ext_X2(1:N-2)),(ext_Y2(3:N) - ext_Y2(1:N-2)));  % centered difference normal, span of 1
    normal_angle(1) = atan2(-(ext_X2(2) - ext_X2(N)),(ext_Y2(2) - ext_Y2(N)));                  % initial point                                  
    normal_angle(N) = atan2(-(ext_X2(1) - ext_X2(N-1)),(ext_Y2(1) - ext_Y2(N-1)));              % final point
    COSINE_TERM = cos(normal_angle);        % direction cosine
    SINE_TERM = sin(normal_angle);          % direction sine

    if N <= minnum  % check that length of ext_X2 is at least minnum before smoothing.
        crossing = 0;
        break
    end


    if i == 1                               % FIRST ITERATION: Analysis of the boundary and increment direction must be chosen.
        ext_X2R = round(ext_X2);            % round points for use as indices
        ext_Y2R = round(ext_Y2);
        index = sub2ind(size(I),ext_Y2R, ext_X2R);  % convert row (Y) and column (X) indices into a linear index
        rind1  = I(index);                  % rind intensity values of the original boundary.
        avgrindref = mean(rind1);
        medrindref = median(single(rind1));
        percrindref = prctile(single(rind1),[0 5 25 50 75 95 100]);


        d = 1;                              % increment to step    
        ext_X2 = d*COSINE_TERM + ext_X2;    % new contour is the old contour plus d*COS or d*SIN
        ext_Y2 = d*SINE_TERM + ext_Y2;
        ext_X2R = round(ext_X2);                % round points for use as indices
        ext_Y2R = round(ext_Y2);
        index = sub2ind(size(I),ext_Y2R, ext_X2R);  % convert row (Y) and column (X) indices into a linear index
        rind2  = I(index);                       % rind intensity values of the shifted boundary.

        avgdiff = mean(single(rind2) - single(rind1));
        if avgdiff < 0                          % if the first offset was in the wrong direction....
            d = -2;                             % reverse shift direction and shift by 2 to remove initial offset and move 1  more in the correct direction
            ext_X2 = d*COSINE_TERM + ext_X2;    % new contour is the old contour plus d*COS or d*SIN
            ext_Y2 = d*SINE_TERM + ext_Y2;
            ext_X2R = round(ext_X2);            % round points for use as indices
            ext_Y2R = round(ext_Y2);
            index = sub2ind(size(I),ext_Y2R, ext_X2R);  % convert row (Y) and column (X) indices into a linear index
            rind2  = I(index);                  % rind intensity values of the shifted boundary.
            d = -1;                             % set d to negative 1 for remaining incrementing steps
        end
        medrind(i) = median(single(rind2)); % save the first value of medrind
        continue;                           % Don't do anything else for i = 0;


    else  % ALL OTHER ITERATION NUMBERS
        ext_X2 = d*COSINE_TERM + ext_X2;    % new contour is the old contour plus d*COS or d*SIN
        ext_Y2 = d*SINE_TERM + ext_Y2;
        ext_X2R = round(ext_X2);            % round points for use as indices
        ext_Y2R = round(ext_Y2);
        index = sub2ind(size(I),ext_Y2R, ext_X2R);  % convert row (Y) and column (X) indices into a linear index
        rind2  = I(index);                  % rind intensity values of the shifted boundary.
        medrind(i) = median(single(rind2));
    end
    % ===================  END STEP 4 ==============================================================================


    % ================= STEP 5 - DECIDING WHEN TO STOP ========================================================
    % Decide when to bail out of the iteration process
    if medrind(i) < medrindref              % >>>>>> BREAK OUT WHEN THE MEDIAN VALUE OF THE CURRENT CONTOUR IS LESS THAN THE MEDIAN OF THE REFERENCE CONTOUR
        crossing = 1;                       % set the crossing variable before breaking out
        break                               % break out!
    end
    % ================= STEP 5 ================================================================================

end
% ===================  END MAIN LOOP  ==============================================================================


% ===================  STEP 6 ==============================================================================
% interpolation of the rindthickness based on last two points only
if crossing == 1    
    ya = i - 1;
    yb = i;
    xa = medrind(i-1);
    xb = medrind(i);
    x = medrindref;
    avgrindthickness = ya + ((yb - ya)/(xb-xa))*(x-xa);
else 
    avgrindthickness = i;
end
% ===================  STEP 6 ==============================================================================


% ===================  STEP 7 ==============================================================================
% RECALCULATE THE INTERIOR CURVES BASED ON INTERPOLATED VALUE
d = d*(-1)*(i - avgrindthickness);      % d gives current direction, (-1) reverses that, and (i-avgrind) is how far to travel back);
int_X = d*COSINE_TERM + ext_X2;         % new contour is the old contour plus d*COS or d*SIN
int_Y = d*SINE_TERM + ext_Y2;
int_X = double(int_X);
int_Y = double(int_Y);
% ===================  STEP 7 ==============================================================================


% =================== STEP 8 ==============================================================================
% PLOTTTING (OPTIONAL)
if plotting == 1 
    figure;
    imshow(I)
    hold on
    plot(ext_X,ext_Y,'y')
    plot(int_X,int_Y,'y')
end
% =================== STEP 8 ==============================================================================

end
