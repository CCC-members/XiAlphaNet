% try and erease
load('Data\Model_Parameters\parameters.mat')
tic
[x] =  Xi_ALphaNET_delays(parameters)
toc

%%
%% 1) Generate random data
% Let's generate a 5x5 matrix of random integers between 1 and 50
B = randi(50, 5, 5);

% Pick a random integer x between 1 and 50
x = randi(50);

% Display them
disp('Matrix B =');
disp(B);
fprintf('x = %d\n\n', x);

%% 2) Compute the absolute differences
diffMatrix = abs(B - x);

%% 3) Find the position where |x - B| is minimized
[minVal, linearIdx] = min(diffMatrix(:));  % minVal is the minimum difference
[row, col] = ind2sub(size(B), linearIdx);  % row and col give the 2D location

fprintf('Minimum absolute difference: %d\n', minVal);
fprintf('Occurs at position (row=%d, col=%d)\n', row, col);
fprintf('The value in B at that position is B(%d, %d) = %d\n\n', row, col, B(row, col));

%% 4) Verify correctness of the solution
% The easiest check is to see if there is any difference smaller than minVal.
% We'll compare minVal to the minimum of all differences, which should match.
allDiffs = diffMatrix(:);             % Flatten all differences
checkMin = min(allDiffs);             % Should be identical to minVal
fprintf('Check: The smallest difference among all entries is %d\n', checkMin);

if abs(minVal - checkMin) < 1e-12
    disp('Verification passed: minVal is indeed the smallest absolute difference.');
else
    disp('Verification failed: something went wrong.');
end
