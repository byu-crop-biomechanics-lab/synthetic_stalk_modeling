function writespline(len,data,k,outputName)
    %define empty spline and number of x-y points
    spline = '';

    %run through 1-column arrays of the x and y data points for the spline, and add to the end of the string with the correct formatting
    for i = 1:len 
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), '); 
    end

    %write string to file using the specimen name
    I = sprintf('%d',k);
    fid = fopen(strcat(outputName,I,'.txt'),'wt');
    fprintf(fid, spline);
    fclose(fid);
end