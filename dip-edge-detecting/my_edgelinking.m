function output = my_edgelinking(input, rstart, cstart)
%in this function, you should finish the edge linking utility.
%the input parameters are a matrix of a binary image containing the edge
%information and coordinates of one of the edge points of a obeject
%boundary, you should run this function multiple times to find different
%object boundaries
%the output parameter is a Q-by-2 matrix, where Q is the number of boundary
%pixels. B holds the row and column coordinates of the boundary pixels.
%you can use different methods to complete the edge linking function
%the better the quality of object boundary and the more the object boundaries, you will get higher scores
global EDGEIM;
global ROWS;
global COLS;

global noPoint;
global thereIsAPoint;
% global lastPoint;
EDGEIM = input ~= 0;                     % make sure image is binary.
EDGEIM = bwmorph(EDGEIM,'clean');     % Remove isolated pixels
EDGEIM = bwmorph(EDGEIM,'skel',Inf);  % and make sure edges are thinned. I
% think using 'skel' is better than 'thin'

EDGEIM = double(EDGEIM);   % Convert to double to allow the use of -ve labelings
[ROWS, COLS] = size(EDGEIM);
edgeNo = 1;

noPoint = 0;
thereIsAPoint = 1;
% lastPoint = 2;

edgepoints = [rstart cstart];      % Start a new list for this edge.
EDGEIM(rstart,cstart) = -edgeNo;   % Edge points in the image are
% encoded by -ve of their edgeNo.

[status, r, c] = nextpoint(rstart,cstart, edgeNo); % Find next connected
% edge point.

while status ~= noPoint
    edgepoints = [edgepoints             % Add point to point list
        r    c   ];
    EDGEIM(r,c) = -edgeNo;               % Update edge image
    
    [status, r, c] = nextpoint(r,c, edgeNo); % Otherwise keep going
end

% Now track from original point in the opposite direction

% First reverse order of existing points in the edge list
edgepoints = flipud(edgepoints);

% ...and start adding points in the other direction.
[status, r, c] = nextpoint(rstart,cstart, edgeNo);

while status ~= noPoint
    edgepoints = [edgepoints
        r    c   ];
    EDGEIM(r,c) = -edgeNo;
        [status, r, c] = nextpoint(r,c, edgeNo);
end

% Final check to see if this edgelist should have start and end points
% matched to form a loop.  If the number of points in the list is four or
% more (the minimum number that could form a loop), and the endpoints are
% within a pixel of each other, append a copy if the first point to the
% end to complete the loop

if length(edgepoints) >= 4
    if abs(edgepoints(1,1) - edgepoints(end,1)) <= 1  &&  ...
            abs(edgepoints(1,2) - edgepoints(end,2)) <= 1
        edgepoints = [edgepoints
            edgepoints(1,:)];
    end
end

output = edgepoints;

function [status, nextr, nextc] = nextpoint(rp,cp, edgeNo)

global EDGEIM;
global ROWS;
global COLS;
global noPoint;
global thereIsAPoint;

% row and column offsets for the eight neighbours of a point
roff = [-1  0  1  0 -1 -1  1  1];
coff = [ 0  1  0 -1 -1  1  1 -1];

r = rp+roff;
c = cp+coff;

% Find indices of arrays of r and c that are within the image bounds
ind = find((r>=1 & r<=ROWS) & (c>=1 & c<=COLS));

% If we get here there were no junction points.  Search through neighbours
% and return first connected edge point that itself has less than two
% neighbours connected back to our current edge.  This prevents occasional
% erroneous doubling back onto the wrong segment

checkFlag = 0;
for i = ind
    if EDGEIM(r(i),c(i)) == 1
        n = neighbours(r(i),c(i));
        if sum(n==-edgeNo) < 2
            nextr = r(i);
            nextc = c(i);
            status = thereIsAPoint;
            return;             % break out and return with the data
            
        else                    % Remember this point just in case we
            checkFlag = 1;      % have to use it
            rememberr = r(i);
            rememberc = c(i);
        end
        
    end
end

% If we get here (and 'checkFlag' is true) there was no connected edge point
% that had less than two connections to our current edge, but there was one
% with more.  Use the point we remembered above.
if checkFlag
    nextr = rememberr;
    nextc = rememberc;
    status = thereIsAPoint;
    return;                % Break out
end

% If we get here there was no connecting next point at all.
nextr = 0;
nextc = 0;
status = noPoint;

% Function to get the values of the 8 neighbouring pixels surrounding a point
% of interest.  The values are ordered from the top-left point going
% anti-clockwise around the pixel.
function n = neighbours(rp, cp)

global EDGEIM;
global ROWS;
global COLS;

% row and column offsets for the eight neighbours of a point
roff = [-1  0  1  1  1  0 -1 -1];
coff = [-1 -1 -1  0  1  1  1  0];

r = rp+roff;
c = cp+coff;

% Find indices of arrays of r and c that are within the image bounds
ind = find((r>=1 & r<=ROWS) & (c>=1 & c<=COLS));

n = zeros(1,8);
for i = ind
    n(i) = EDGEIM(r(i),c(i));
end
