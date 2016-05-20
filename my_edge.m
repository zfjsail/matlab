function output = my_edge(input, method, threshold)
%in this function, you should finish the edge detection utility.
%the input parameter is a matrix of a gray image
%the output parameter is a matrix contains the edge index using 0 and 1
%the entries with 1 in the matrix shows that point is on the edge of the
%image
%you can use different methods to complete the edge detection function
%the better the final result and the more methods you have used, you will get higher scores  
if strcmp(method, 'sobel')
    k = [1 2 1; 0 0 0; -1 -2 -1]; %weight
    H = conv2(input,k, 'same');
    V = conv2(input,k','same');
    E = sqrt(H.*H + V.*V);
    output = logical(E > threshold);
else if strcmp(method, 'prewitt')
        k = [1 1 1; 0 0 0; -1 -1 -1];
        H = conv2(input,k, 'same');
        V = conv2(input,k','same');
        E = sqrt(H.*H + V.*V);
        output = logical(E > threshold);
        else if strcmp(method, 'roberts')
                k = [0 0 0; 0 1 0; 0 0 -1];
                H = conv2(input,k, 'same');
                V = conv2(input,k','same');
                E = abs(H) + abs(V); 
                output = logical(E > threshold);     
            else if strcmp(method, 'canny')
                    T_Low = 0.075;
                    T_High = 0.175;
                    % X-axis direction edge detection
                    filterx=d2dgauss(5, 1, 5, 1, pi/2);
                    Ix= conv2(input,filterx,'same');
                    
                    % Y-axis direction edge detection
                    filtery=d2dgauss(5,1,5,1,0);
                    Iy=conv2(input,filtery,'same');
                    
                    % Norm of the gradient (Combining the X and Y directional derivatives)
                    NVI=sqrt(Ix.*Ix+Iy.*Iy);
                    
                    % Filter for horizontal and vertical direction
                    % this filter has better performance than that in
                    % slides
                    KGx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
                    KGy = [1, 2, 1; 0, 0, 0; -1, -2, -1];
                    %Convolution by image by horizontal and vertical filter
                    Filtered_X = conv2(NVI, KGx, 'same');
                    Filtered_Y = conv2(NVI, KGy, 'same');
                    
                    %Calculate directions/orientations
                    arah = atan2 (Filtered_Y, Filtered_X);
                    arah = arah*180/pi;
                    
                    pan=size(NVI,1);
                    leb=size(NVI,2);
                    
                    %Adjustment for negative directions, making all directions positive
                    for i=1:pan
                        for j=1:leb
                            if (arah(i,j)<0)
                                arah(i,j)=360+arah(i,j);
                            end;
                        end;
                    end;
                    
                    arah2=zeros(pan, leb);
                    
                    %Adjusting directions to nearest 0, 45, 90, or 135 degree
                    for i = 1  : pan
                        for j = 1 : leb
                            if ((arah(i, j) >= 0 ) && (arah(i, j) < 22.5) || (arah(i, j) >= 157.5) && (arah(i, j) < 202.5) || (arah(i, j) >= 337.5) && (arah(i, j) <= 360))
                                arah2(i, j) = 0;
                            elseif ((arah(i, j) >= 22.5) && (arah(i, j) < 67.5) || (arah(i, j) >= 202.5) && (arah(i, j) < 247.5))
                                arah2(i, j) = 45;
                            elseif ((arah(i, j) >= 67.5 && arah(i, j) < 112.5) || (arah(i, j) >= 247.5 && arah(i, j) < 292.5))
                                arah2(i, j) = 90;
                            elseif ((arah(i, j) >= 112.5 && arah(i, j) < 157.5) || (arah(i, j) >= 292.5 && arah(i, j) < 337.5))
                                arah2(i, j) = 135;
                            end;
                        end;
                    end;
                    %Calculate magnitude
                    magnitude = (Filtered_X.^2) + (Filtered_Y.^2);
                    magnitude2 = sqrt(magnitude);
                    
                    BW = zeros (pan, leb);
                    
                    %Non-Maximum Supression
                    for i=2:pan-1
                        for j=2:leb-1
                            if (arah2(i,j)==0)
                                BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i,j+1), magnitude2(i,j-1)]));
                            elseif (arah2(i,j)==45)
                                BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j-1), magnitude2(i-1,j+1)]));
                            elseif (arah2(i,j)==90)
                                BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j), magnitude2(i-1,j)]));
                            elseif (arah2(i,j)==135)
                                BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j+1), magnitude2(i-1,j-1)]));
                            end;
                        end;
                    end;
                    BW = BW.*magnitude2;
                    
                    %Hysteresis thresholding
                    T_Low = T_Low * max(max(BW));
                    T_High = T_High * max(max(BW));
                    
                    T_res = zeros (pan, leb);
                    
                    for i = 1  : pan
                        for j = 1 : leb
                            if (BW(i, j) < T_Low)
                                T_res(i, j) = 0;
                            elseif (BW(i, j) > T_High)
                                T_res(i, j) = 1;
                                %Using 8-connected components
                            elseif ( BW(i+1,j)>T_High || BW(i-1,j)>T_High || BW(i,j+1)>T_High || BW(i,j-1)>T_High || BW(i-1, j-1)>T_High || BW(i-1, j+1)>T_High || BW(i+1, j+1)>T_High || BW(i+1, j-1)>T_High)
                                T_res(i,j) = 1;
                            end;
                        end;
                    end;
                    output = logical(T_res);
                else if strcmp(method, 'log')
                        if ~exist('threshold', 'var') || isempty(threshold)
                            threshold = 0;
                        end

                            % X-axis direction edge detection

                            % processing individually is better than
                            % combine step 1 & 2 together
                            filterx=d2dgauss(5, 1, 5, 1, pi/2);
                            Ix= conv2(input,filterx,'same');
                            
                            % Y-axis direction edge detection
                            filtery=d2dgauss(5,1,5,1,0);
                            Iy=conv2(input,filtery,'same');
                            
                            % Norm of the gradient (Combining the X and Y directional derivatives)
                            NVI=sqrt(Ix.*Ix+Iy.*Iy);
                            
                            mask = [1 1 1; 1 -8 1; 1 1 1];
                            O = conv2(NVI, mask, 'same');
                            
                        output = zeroCross(O, threshold);
                    end
                end
            end
        end
    end
end

function h = d2dgauss(n1, std1, n2, std2, theta)
r = [cos(theta)  -sin(theta); sin(theta)  cos(theta)];
h = zeros(n2, n1);
for i = 1 : n2
    for j = 1 : n1
        u = r * [j - (n1 + 1)/2 i-(n2 + 1)/2]';
        h(i, j) = gauss(u(1), std1)*gauss(u(2), std2);
    end
end
h = h/sqrt(sum(sum(h.*h)));
end

function y = gauss(x, std)
y = exp(-x^2/(2*std^2))/(std*sqrt(2*pi));
end

function idx = zeroCross(input, threshold, direction)
if ~ismatrix(input)
    error('The input must be a two dimensional array.');
end
if ~exist('threshold', 'var') || isempty(threshold)
    threshold = 0;
end
if ~exist('direction', 'var') || isempty(direction)
    % all directions
    idx = zeroCross(input, threshold, 'horizontal');
    idx = idx | zeroCross(input, threshold, 'vertical');
    idx = idx | zeroCross(input, threshold, '45');
    idx = idx | zeroCross(input, threshold, '135');
    return;
end
if strcmp(direction, 'horizontal')
    mask = [0  0  0;
            -1  0  1;
            0  0  0];
elseif strcmp(direction, 'vertical')
    mask = [0  -1  0;
            0  0  0 ;
            0  1  0];
elseif strcmp(direction, '135')
    mask = [-1  0 0;
            0  0 0;
            0  0 1];
elseif strcmp(direction, '45')
    mask = [0  0  1;
            0  0  0;
           -1  0  0];
end
%% check if there is a change in sign
s = sign(input);
t = conv2(s, mask, 'same');
idx = (abs(t) == 2);
% To consider 0 for changes in sign, use (abs(t) > 0)

%% thresholding
if threshold > 0
    a = conv2(input, mask, 'same');
    idx = idx & (abs(a) > threshold);
end
end