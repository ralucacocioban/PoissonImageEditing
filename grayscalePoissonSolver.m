function [result] = grayscalePoissonSolver(source, mask)

[srcW srcH channels] = size(source);    

if(channels ~= 1)
    fprintf('Image is not grayscale. \n');
end
% --------------------------------------------
% get the nb of pixels that are within the mask
% --------------------------------------------

n = sum(mask(:) == 1);

indexes = getIndexes(mask, srcH, srcW);

if(sum(mask(:) == 1) ~= n)
    fprintf('Index and coefficient matrix size mismatch. \n');
end

%---------------------------------------------
% build the coefficient matrix and the solution vector
%---------------------------------------------

coefficientMatrix = getCoefMatrix(srcH, srcW, mask, indexes);
solVector = getSolutionVect(srcH, srcW, source, mask);

%---------------------------------------------
% solve for the sparse matrix
%---------------------------------------------
solution = coefficientMatrix\( solVector(1,:)' );


result = reconstructImg(solution, mask, source, indexes);

end

function resultImg = reconstructImg(solution, mask, source, indexes)

[sourceW sourceH channels] = size(source);

resultImg = source;

for x = 1:sourceW
    for y = 1:sourceH
        if mask(x, y) ~= 0
            index = indexes(x, y);
            resultImg(x, y) = solution(index);
        end
    end
end
    
end

function coefM = getCoefMatrix(srcH, srcW, mask, indexes)

countInsidePixels = 0;

n = sum(mask(:) == 1);
%---------------------------------------------
% coefficient matrix
%---------------------------------------------
coefM = spalloc(n, n, 3*n);

% go only through the elements inside the mask
% for neighbor pixels within the mask the coeff value = -1
% for neighbor pixels outside the mask coeff value = 0
% Np otherwise


for i = 2:srcW-1
    for j = 2:srcH-1

        if mask(i, j) ~= 0
            countInsidePixels = countInsidePixels + 1;   
            
            % if on the top
            if mask(i-1, j) ~= 0 
                coefM(countInsidePixels, indexes(i-1, j)) = -1;
            end
            
            % if on the left
            if mask(i, j-1) ~= 0
                coefM(countInsidePixels, indexes(i, j-1)) = -1;
            end
            
            % if on the bottom            
            if mask(i+1, j) ~= 0
                coefM(countInsidePixels, indexes(i+1, j)) = -1;
            end
            
            % if on the right side
            if mask(i, j+1) ~= 0
                coefM(countInsidePixels, indexes(i, j+1)) = -1;
            end       
            
            coefM(countInsidePixels, countInsidePixels) = 4;
        end
    end
end

end


function indexes = getIndexes(mask, sourceH, sourceW)

% indicator of row index in the solution vector
indexes = zeros(sourceW, sourceH);

insidePixels = 0;

for i = 1:sourceW
    for j = 1:sourceH
        if (mask(i,j) ~= 0)
            insidePixels = insidePixels + 1;            
            indexes(i, j) = insidePixels;
        end
    end
end

end

function solVector = getSolutionVect(srcH, srcW, source, mask)
    
  n = sum(mask(:) == 1);
  solVector = zeros(1, n);

  countInsidePixels = 0;
  
  for i = 2:srcW-1
    for j = 2:srcH-1

        if mask(i, j) ~= 0
            countInsidePixels = countInsidePixels + 1;   
            
            
            if mask(i-1, j) == 0 
                solVector(1, countInsidePixels) = solVector(1, countInsidePixels) + source(i-1, j);
            end
            
            
            if mask(i, j-1) == 0
                solVector(1, countInsidePixels) = solVector(1, countInsidePixels) + source(i, j-1);
            end
            
                       
            if mask(i+1, j) == 0
                solVector(1, countInsidePixels) = solVector(1, countInsidePixels) + source(i+1, j);
            end
            
            if mask(i, j+1) == 0
                solVector(1, countInsidePixels) = solVector(1, countInsidePixels) + source(i, j+1);
            end       
        end
    end
  end
end
