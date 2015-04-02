img =  imread('sun.jpg') ;
source = double(img);
gray_source = rgb2gray(img); 

mask = roipoly(img);

destination = double(imread('bungalow.jpg'));

[srcH srcW] = size(source);
[destH destW] = size(destination);
[maskH maskW] = size(mask);

% check if the mask fits into the destination image
if (maskH > destH || maskW > destW || maskH > srcH || maskW > srcW)
    fprintf('The mask exceeds the destination size \n');  
else   
    
    
    result = seamlessMixedCloningPoisson(source, destination, mask);
    figure('Name', 'Result of mixed Poisson seamless cloning');
    imshow(mat2gray(result));
    
    resultSimple = seamlessCloningPoisson(source, destination, mask);
    figure('Name', 'Result of simple Poisson seamless cloning');
    imshow(mat2gray(resultSimple));
   
    %{
    
    task2Result = simplePoissonSolver(source, mask);
    figure('Name', 'Result of simple Poisson solver for RGB images');
    imshow(mat2gray(task2Result));
    
    task2GraycsaleResult = grayscalePoissonSolver(gray_source, mask);
    figure('Name', 'Result of simple Poisson solver for grayscale images');
    imshow(mat2gray(task2GraycsaleResult));
    
    %}
    
    
end
