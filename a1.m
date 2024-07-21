function set12 = a1_20292366(elist) 
% Replace 00000000 with your student number and rename file
% Function for CISC271, Winter 2022, Assignment #1
%
% IN:
%     elist - Mx2 array of edges, each row is a pair of vertices
% OUT:
%     set12 - Nx1 vertex clustering, -1 for SET1 and +1 for SET2

    % Problem size: number of vertices in the graph
    n = max(elist(:));

    % %
    % % STUDENT CODE GOES HERE: replace this trivial adjacency matrix
    % %     with your computation; the first line is a hint for how to
    % %     initialize your matrix correctly
    % %
    
    %   This initializes a matrix of size n, where n is the maximum no. of
    %   vertices in the edge set.
    
    A = zeros(n);

    %   This for-loop is used to create an adjacency matrix using the given
    %   edge set stored in elist as an nx2 array
    %   It does so by traversing all the elements stored in the edge list, 
    %   and replacing the zeroes in the cells corresponding to vertices
    %   connected by an edge in the initialized matrix with ones
    
    for i = 1:length(elist(:,1))
        A(elist(i,1),elist(i,2)) = 1;
    end

    % %
    % % STUDENT CODE GOES HERE: replace this constant clustering vector
    % %     with your computation that uses the Fiedler vector; the
    % %     vector SET12 should be plus/minus 1 for automated grading
    % %    

    %   Since we know from A1 tutorial.mlx that the "degree matrix is a 
    %   diagonal matrix where entry (i,i) is the number of edges connected 
    %   to vertex i", we create a column vector that is the sum of each row
    %   and then change it to be along the diagonal of a 20x20 matrix

    D = diag(sum(A,2));
    L = D - A;

    %   From A1 Tutorial.mlx - Recall: Laplacian matrix is 
    %   (degree matrix) - (adjacency matrix) - We find this to calculate
    %   the eigenvectors for the graph
    
    [eigVec, ~] = eig(L);

    %   Here we calculate the Fiedler Vector using the second column of the
    %   eigenvector matrix which helps us create a clustering that we
    %   further modify below

    fiedlerVec = eigVec(:,2);

    %   The set12 is a combination of 2 sets - one that has a value of +1
    %   and one that equals zero, and since zero is neither strictly
    %   positive or negative, we transform all the zeroes to -1, giving us
    %   a set12 with +1 and -1 values as our final combined set.
    
    set12 = fiedlerVec <=0;
    set12 = (set12*2) - 1;

    % %
    % % STUDENT CODE GOES HERE: replace this trivial display to console
    % %     with your computation from the vector SET12
    % %
    
    %   We initialize sets 1 and 2 where we want to separate the indices of
    %   values -1 and 1 respectively. We know each cluster is going to have
    %   10 vertices since this is a sort of binary clustering which is only
    %   technically effective for 2 dimensional data

    set1 = zeros(1,10);
    set2 = zeros(1,10);
    
    %   A for-loop is used to traverse through set12 and extract the indices
    %   of values with -1 to set 1 and +1 to set 2.
    
    for i = 1:length(set12)
        if set12(i) == -1, set1(i) = i;
        elseif set12(i) == 1, set2(i) = i;
        end
    end

    %   Here, all the stray zero values left over from the initialization
    %   of sets 1 and 2 are discarded to create a clear output

    set1(set1 == 0) = [];
    set2(set2 == 0) = [];

    %   An optimized version of this code can be found below where time
    %   complexity is reduced since there isnt any need for the loop and
    %   zeroes are automatically erased.
    %   I wasn't sure whether this would be marked on the basis of there
    %   being a for-loop or not since it is one of the learning outcomes on
    %   the assignment -

    %   set1 = find(set12 == -1);
    %   set2 = find(set12 == 1);
    

    %   This prints out the 2 clusters of vertices as clear output

    disp('Set 1 vertices are:');
    disp(set1);
    disp('Set 2 vertices are:');
    disp(set2);

    % Plot the graph, Cartesian and clustered
    plot271a1(A, set12);
end

function plot271a1(Amat, cvec)
% PLOTCLUSTER(AMAT,CVEC) plots the adjacency matrix AMAT twice;
% first, as a Cartesian grid, and seconnd, by using binary clusters
% in CVEC to plot the graph of AMAT based on two circles
%
% INPUTS: 
%         Amat - NxN adjacency matrix, symmetric with binary entries
%         cvec - Nx1 vector of class labels, having 2 distinct values
% OUTPUTS:
%         none
% SIDE EFFECTS:
%         Plots into the current figure

    % %
    % % Part 1 of 2: plot the graph as a rectangle
    % %

    % Problem size
    [m n] = size(Amat);

    % Factor the size into primes and use the largest as the X size
    nfact = factor(n);
    nx = nfact(end);
    ny = round(n/nx);

    % Create a grid and pull apart into coordinates; offset Y by +2
    [gx, gy] = meshgrid((1:nx) - round(nx/2), (1:ny) + 2);

    % Offset the odd rows to diagram the connections a little better
    for ix=1:2:ny
        gx(ix, :) = gx(ix, :) + 0.25*ix;
    end

    % The plot function needs simple vectors to create the graph
    x = gx(:);
    y = flipud(gy(:));

    % Plot the graph of A using the Cartesian grid
    plot(graph(tril(Amat, -1), 'lower'), 'XData', x, 'YData', y);
    axis('equal');

    % %
    % % Part 2 of 2: plot the graph as pair of circles
    % %
    % Set up the X and Y coordinates of each graph vertex
    xy = zeros(2, numel(cvec));

    % Number of cluster to process
    kset = unique(cvec);
    nk = numel(kset);

    % Base circle is radius 2, points are centers of clusters
    bxy = 2*circlen(nk);

    % Process each cluster
    for ix = 1:nk
        jx = cvec==kset(ix);
        ni = sum(jx);
        xy(:, jx) = bxy(:, ix) + circlen(ni);
    end

    hold on;
    plot(graph(Amat), 'XData', xy(1,:), 'YData', xy(2,:));
    hold off;
    title(sprintf('Clusters of (%d,%d) nodes', ...
        sum(cvec==kset(1)), sum(cvec==kset(2))));
end

function xy = circlen(n)
% XY=CIRCLEN(N) finds N 2D points on a unit circle
%
% INPUTS:
%         N  - positive integer, number of points
% OUTPUTS:
%         XY - 2xN array, each column is a 2D point

    xy = [cos(2*pi*(0:(n-1))/n) ; sin(2*pi*(0:(n-1))/n)];
end
