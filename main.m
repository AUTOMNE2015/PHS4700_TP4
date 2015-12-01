function main
     %little matlab fairy with great magic

    % simulation 1
    positionInit = [-10 -10 15];
    option = 1;

%     % simulation 2
%     positionInit = [13 10 25];
%     option = 1;

%     % simulation 3
%     positionInit = [-10 -10 15];
%     option = 2;
    
%     % simulation 4
%     positionInit = [13 10 25];
%     option = 2;

    % generer les directions
    nbreDePoints = 30;
    coinBloc = [-2 -2 1];    
    directions = zeros(3);
    count = 1;
    incrementX = 11 / nbreDePoints;
    incrementY = 11 / nbreDePoints;
    incrementZ = 23 / (nbreDePoints*2);
    
    for i = 1:nbreDePoints
        for j = 1:nbreDePoints
            for k = 1:(2*nbreDePoints)
                directions(count, :) = [coinBloc(1)+i*incrementX coinBloc(2)+j*incrementY  coinBloc(3)+k*incrementZ]- positionInit;
                directions(count, :) = directions(count, :) / norm(directions(count, :));
                count = count + 1;
            end
        end
    end
    tracerPoints(positionInit, directions, nbreDePoints, option);   
end

%varaibel global
function setGlobalpremier(val)
global premier
premier = val;
end

function r = getGlobalpremier
global premier
r = premier;
end 

function drawBigCube
     A = [0 0 5];%[0 0 0];
     B = [7 0 5];%[1 0 0];
     C = [0 7 5];%[0 1 0];
     D = [0 0 20];%[0 0 1];
     E = [0 7 20];%[0 1 1];
     F = [7 0 20];%[1 0 1];
     G = [7 7 5];%[1 1 0];
     H = [7 7 20];%[1 1 1];
     P = [A;B;F;H;G;C;A;D;E;H;F;D;E;C;G;B];
     h = plot3(P(:,1),P(:,2),P(:,3))
     h.Color = 'black';
     
     A = [3 3 12];%[0 0 0];
     B = [4 3 12];%[1 0 0];
     C = [3 5 12];%[0 1 0];
     D = [3 3 17];%[0 0 1];
     E = [3 5 17];%[0 1 1];
     F = [4 3 17];%[1 0 1];
     G = [4 5 12];%[1 1 0];
     H = [4 5 17];%[1 1 1];
     P = [A;B;F;H;G;C;A;D;E;H;F;D;E;C;G;B];
     h = plot3(P(:,1),P(:,2),P(:,3))
     h.Color = 'black';
     axis([-30, 30, -30, 30, -30, 30]);
end




%choisir le secteur angulaire de départ (N directions de lignes)
%tracer les lignes
    %verifier si elle pénètre ou si elle est réfléchie par la surface
        %calculer refraction
        %calculer reflexion
%identifier les rayons qui atteignent l'objet
%tracer la position de l'image virtuelle

function y = delta()
    y = 0.01;
end

function y = positionBloc()
    y =[0 0 20;
        0 7 20;
        7 7 20;
        7 0 20;
        0 0 5;
        0 7 5;
        7 7 5;
        7 0 5];
end

function y = normaleSurface()
    y = [-1 0 0;
        0 1 0;
        1 0 0;
        0 -1 0;
        0 0 1;
        0 0 -1];
end

function y = positionPetitBloc()
    y = [3 3 17;
        3 5 17;
        4 5 17;
        4 3 17;
        3 3 12;
        3 5 12;
        4 5 12;
        4 3 12];
end

function y = positionGrosseBoite()
    y = [-30 -30 30;
        -30 30 30;
        30 30 30;
        30 -30 30;
        -30 -30 -30;
        -30 30 -30;
        30 30 -30;
        30 -30 -30];
end

function y = indiceRefraction(position, direction, option)
    posBloc = positionBloc();
    
%     for i = 1:6
%         planCourant = posBloc(i);  
%         vecteurNormal = cross(planCourant(1),planCourant(2));
%         vecteurNormalUnitaire = vecteurNormal/norm(vecteurNormal);
%     end
    
    % arriere
    posArriere = [position(1) - delta()*direction(1)
                  position(2) - delta()*direction(2)
                  position(3) - delta()*direction(3)];  
    % devant
    posDevant = [position(1) + delta()*direction(1)
                 position(2) + delta()*direction(2)
                 position(3) + delta()*direction(3)];
             
    min = posBloc(5, :);
    max = posBloc(3, :);
    %n1
    if(option == 1)
        n1 = 1.0;
    else
        n1 = 1.33;
    end
    if(min(1) < posArriere(1) && max(1) > posArriere(1))
        if(min(2) < posArriere(2) && max(2) > posArriere(2))
             if(min(3) < posArriere(3) && max(3) > posArriere(3))
                 %interieur
                 if(option == 1)
                     n1 = 1.5;
                 else
                     n1 = 1.1;
                 end
             end
        end
    end
    %n2
    if(option == 1)
        n2 = 1.0;
    else
        n2= 1.33;
    end
    if(min(1) < posDevant(1) && max(1) > posDevant(1))
        if(min(2) < posDevant(2) && max(2) > posDevant(2))
             if(min(3) < posDevant(3) && max(3) > posDevant(3))
                 %interieur
                 if(option == 1)
                     n2 = 1.5;
                 else
                     n2 = 1.1;
                 end
             end
        end
    end
    y = [n1 n2];
end

function y = isBetweenTwoPoints(point1, point2, point3)
    y =0;
    if(numel(point1)==0)
        return;
    end
    max = point3;
    min = point2;
    for i = 1:3
        if(point2(i) > point3(i))
            max(i) = point2(i);
            min(i) = point3(i);
        end
    end
    
    if(point1(1) <= max(1) && point1(1) >= min(1))
        if(point1(2) <= max(2) && point1(2) >= min(2))
             if(point1(3) <= max(3) && point1(3) >= min(3))
                 y = 1;
             end
        end
    end          
end


function tracerPoints(position, directions, nombrePoints, option) 
    %directions = tableau de vecteur
    %N fois tracerUneLigne
    hold on
    scatter3(position(1),position(2),position(3));
    %drawBigCube();
    total = nombrePoints*nombrePoints*2*nombrePoints;
    for i = 1:total
        clc();
        fprintf('Drawing...\n%5.2f %%\n', (i/total)*100);
        
        setGlobalpremier(1); %premiere iteration pour le prochain rayon
        detailPoint = tracerUneLigne(position, directions(i, :), option);
        if(detailPoint(1) > 1)
            point = detailPoint(2)*directions(i, :) + position;
            scatter3(point(1),point(2),point(3),10,getColor(detailPoint(1)),'filled'); %TODO : add color LOL
            %fprintf('color %i', detailPoint(1));
        end
    end
    axis([-30, 30, -30, 30, -30, 30]);
end

function y = getColor(noSurface)

    switch(noSurface)
        case 2
            y = [1 0 0];
        case 3
            y = [1 1 0];
        case 4
            y = [0 1 1];
        case 5
            y = [0 1 0];
        case 6
            y = [1 0 1];
        case 7
            y = [0 0 1];
        otherwise
            y = [0 0 0];
    end
%     switch(noSurface)
%         case 2
%             y = 'red';
%         case 3
%             y = 'yellow';
%         case 4
%             y = 'cyan';
%         case 5
%             y = 'green';
%         case 6
%             y = 'pink';
%         case 7
%             y = 'blue';
%         otherwise
%             y = 'black';
%       end
end
%par de la position dans la direction
%calcule le point de collision
%calcule la distance
%tcheck type de collision
    %externe (fini)
    %objet (fini + couleur)
    %bloc 
        %calculer nouvelle direction
        %recursif


function y = tracerUneLigne(position, direction, option)
    
    %calculer collision
    collision = calculerCollision(position, direction);
    typeCollision = collision(1);
    positionCollision = collision(2:4);
    normal = collision(5:7);
    estPremier = getGlobalpremier();
    if(estPremier == 0)
        normal = -normal;  
    end
    
    %calculer distance
    distance = calculerdistance(position, positionCollision);
    
    if(typeCollision == 1)
        angle = calculerAngle2Vecteur(direction, normal);

        n = indiceRefraction(positionCollision, direction, option);
        critique = verifierAngleCritique(angle, n(1), n(2));
        %critique = 1 = réflexion
        %critique = 0 = refraction
        if (critique == 0)
            %fprintf('nouvelle refraction\n');
            nouvelleDirection = calculerDirectionRefraction(direction, normal, n(1), n(2));
        elseif (critique == 1)
           % fprintf('nouvelle reflexion\n');
            nouvelleDirection = calculerDirectionReflexion(direction, normal);
        else
            %ne devrait pas arriver
            y = [0 distance];
            return;
        end
        %on ne sera plus jamais a la premiere iteration pour ce rayon
        setGlobalpremier(0);
        pointCollision = tracerUneLigne(positionCollision, nouvelleDirection, option);
    else
        %fprintf('fin du trace\n');
        pointCollision = [typeCollision distance];
    end
    
    y = [pointCollision(1) (pointCollision(2)+distance)];
end

function y = calculerdistance(position1, position2)
    vecteurDistance = position2 - position1;
    distance = sqrt(vecteurDistance(1)*vecteurDistance(1) + vecteurDistance(2)*vecteurDistance(2) + vecteurDistance(3)*vecteurDistance(3));
    y = distance;
end

function y = calculerCollision(position, direction)
    
    Bloc = positionBloc();
    petitBloc = positionPetitBloc();
    vraiIntersection = [0 0 0];
    minDistance = 1000000000;
    normalSurface = [0 0 0];
    typeCollision = 0;
    % check des collisions avec les face laterales des deux blocs
    for i = 1:4
        
        intersection = intersectLinePlane([position direction], [Bloc(i,:) Bloc(i+1, :) Bloc(i+4, :)]);

        %intersection avec le bloc de verre (ou de polymere bizzare...)
        if(isBetweenTwoPoints(intersection, Bloc(mod(i,4) + 1, :), Bloc(i+4,:)) == 1)
            temp = intersection - position;
            j=1;
            while(temp(j) ==0 && j < 3)
                j = j+1;
            end

            a = temp(j)/direction(j);
            if(a > 0 && norm(temp) < minDistance)
                minDistance = norm(temp);
                vraiIntersection = intersection;
                normalTemp = normaleSurface();
                normalSurface = normalTemp(i, :);
                typeCollision = 1;

                %x
                % garder la surface la plus proche (plus petit x positif)
            end
        end

        intersectionPetit = intersectLinePlane([position direction], [petitBloc(i,:) petitBloc(i+1, :) petitBloc(i+4, :)]);
        %intersection avec le bloc d'acier
        if(isBetweenTwoPoints(intersectionPetit, petitBloc(mod(i,4) + 1, :), petitBloc(i+4,:)) == 1)
            temp = intersectionPetit - position;
            j=1;
            while(temp(j) ==0 && j < 3)
                j = j+1;
            end

            a = temp(j)/direction(j);
            if(a > 0 && norm(temp) < minDistance)
                minDistance = norm(temp);
                vraiIntersection = intersectionPetit;
                normalTemp = normaleSurface();
                normalSurface = normalTemp(i, :);
                typeCollision = i + 1;

                %x
                % garder la surface la plus proche (plus petit x positif)
            end
        end
    end

    % intersection pour les faces superieurs
    intersection = intersectLinePlane([position direction], [Bloc(1,:) Bloc(2, :) Bloc(3, :)]);

    %intersection avec le bloc de verre (ou de polymere bizzare...)
    if(isBetweenTwoPoints(intersection, Bloc(1, :), Bloc(3,:)) == 1)
        temp = intersection - position;
        j=1;
        while(temp(j) ==0 && j < 3)
            j = j+1;
        end

        a = temp(j)/direction(j);
        if(a > 0 && norm(temp) < minDistance)
            minDistance = norm(temp);
            vraiIntersection = intersection;
            normalTemp = normaleSurface();
            normalSurface = normalTemp(5, :);
            typeCollision = 1;
            %x
            % garder la surface la plus proche (plus petit x positif)
        end
    end

    intersectionPetit = intersectLinePlane([position direction], [petitBloc(1,:) petitBloc(2, :) petitBloc(3, :)]);
    %intersection avec le bloc d'acier
    if(isBetweenTwoPoints(intersectionPetit, petitBloc(1, :), petitBloc(3,:)) == 1)
        temp = intersectionPetit - position;
        j=1;
        while(temp(j) ==0 && j < 3)
            j = j+1;
        end

        a = temp(j)/direction(j);
        if(a > 0 && norm(temp) < minDistance)
            minDistance = norm(temp);
            vraiIntersection = intersectionPetit;
            normalTemp = normaleSurface();
            normalSurface = normalTemp(5, :);
            typeCollision = 6;
            %x
            % garder la surface la plus proche (plus petit x positif)
        end
    end


    % intersection pour les faces inferieurs
    intersection = intersectLinePlane([position direction], [Bloc(5,:) Bloc(6, :) Bloc(7, :)]);

    %intersection avec le bloc de verre (ou de polymere bizzare...)
    if(isBetweenTwoPoints(intersection, Bloc(5, :), Bloc(7,:)) == 1)
        temp = intersection - position;
        j=1;
        while(temp(j) ==0 && j < 3)
            j = j+1;
        end

        a = temp(j)/direction(j);
        if(a > 0 && norm(temp) < minDistance)
            minDistance = norm(temp);
            vraiIntersection = intersection;
            normalTemp = normaleSurface();
            normalSurface = normalTemp(6, :);
            typeCollision = 1;
            %x
            % garder la surface la plus proche (plus petit x positif)
        end
    end

    intersectionPetit = intersectLinePlane([position direction], [petitBloc(5,:) petitBloc(6, :) petitBloc(7, :)]);
    %intersection avec le bloc d'acier
    if(isBetweenTwoPoints(intersectionPetit, petitBloc(1, :), petitBloc(3,:)) == 1)
        temp = intersectionPetit - position;
        j=1;
        while(temp(j) ==0 && j < 3)
            j = j+1;
        end

        a = temp(j)/direction(j);
        if(a > 0 && norm(temp) < minDistance)
            minDistance = norm(temp);
            vraiIntersection = intersectionPetit;
            normalTemp = normaleSurface();
            normalSurface = normalTemp(5, :);
            typeCollision = 7;
            %x
            % garder la surface la plus proche (plus petit x positif)
        end
    end

    % refaire la verification avec 2 faces restantes
    %trouver le plus petit x > 0
    % retourner la surface
    positionCollision = vraiIntersection;
    y = [typeCollision, positionCollision, normalSurface];
end

function y = calculerAngle2Vecteur(Vecteur1, Vecteur2)
    CosTheta = dot(Vecteur1,Vecteur2)/(norm(Vecteur1)*norm(Vecteur2));
    y = acos(CosTheta);
end

function y = calculerDirectionRefraction(direction, normale, n1, n2)
    % voir slide p48
    j = cross(direction, normale)/norm(cross(direction, normale));
    k = cross(normale, j)/norm(cross(normale, j)); 
    s = (n1/n2)*(dot(direction,k)); 
    dir = (-normale * sqrt(1-(s^2))) + k*s;
    y = dir/norm(dir); %unitaire
end

function y = calculerDirectionReflexion(direction, normale)
    %voir slide p.35
    res = direction - 2*normale*dot(direction, normale); 
    y = res/norm(res);
end

function y = verifierAngleCritique(angle, n1, n2)
    %formule AngleCritique
    sb = sin(angle) * n1 / n2;
    if(sb > 1)
        y = 1;
    elseif(sb < 1)
        y = 0; 
    else
        y =-1;
        fprintf('paraleleleleleel');
    end
    %-1 = parallele
    %1 = réflexion
    %0 = refraction
end

function point = intersectLinePlane(line, plane)
%trouver vecteur correspondant au plan a laide de ses points
temp1 = plane(:, 4:6)-plane(:,1:3);
temp2 = plane(:, 7:9) - plane(:,1:3);

plane(:,4:6) = temp1;
plane(:,7:9) = temp2;

% n = vectorCross3d(plane(:,4:6), plane(:,7:9));
% 
% if(dot(n,line(:,4:6)) ~= 0)
%     syms d;
%     equationD = n(1)*plane(:,1) + n(2)*plane(:,2) + n(3)*plane(:,3) + d == 0;
%     solutionD = solve(equationD,d);
%     
%     syms t;
%     equationT = n(1)*(line(:,1) + line(:,4)*t) + n(2)*(line(:,2) + line(:,5)*t) + n(3)*(line(:,3) + line(:,6)*t) + solutionD == 0;
%     solutionT = solve(equationT, t);
%     
%     point = [0 0 0];
%     point(1) = line(:,1) + line(:,4)*solutionT;
%     point(2) = line(:,2) + line(:,5)*solutionT;
%     point(3) = line(:,3) + line(:,6)*solutionT;
% else
%     point = [];
% end
   tol = 1e-14;

% unify sizes of data
nLines  = size(line, 1);
nPlanes = size(plane, 1);

% N planes and M lines not allowed 
if nLines ~= nPlanes && min(nLines, nPlanes) > 1
    error('MatGeom:geom3d:intersectLinePlane', ...
        'Input must have same number of rows, or one must be 1');
end

% plane normal
n = vectorCross3d(plane(:,4:6), plane(:,7:9));

% difference between origins of plane and line
dp = bsxfun(@minus, plane(:, 1:3), line(:, 1:3));

% dot product of line direction with plane normal
denom = sum(bsxfun(@times, n, line(:,4:6)), 2);

% relative position of intersection point on line (can be inf in case of a
% line parallel to the plane)
t = sum(bsxfun(@times, n, dp),2) ./ denom;

% compute coord of intersection point
point = bsxfun(@plus, line(:,1:3),  bsxfun(@times, [t t t], line(:,4:6)));

point = round(point,3);
% set indices of line and plane which are parallel to NaN
par = abs(denom) < tol;
point(par,:) = NaN;
end

function c = vectorCross3d(a,b)
%VECTORCROSS3D Vector cross product faster than inbuilt MATLAB cross.
%
%   C = vectorCross3d(A, B) 
%   returns the cross product of the 3D vectors A and B, that is: 
%       C = A x B
%   A and B must be N-by-3 element vectors. If either A or B is a 1-by-3
%   row vector, the result C will have the size of the other input and will
%   be the  concatenation of each row's cross product. 
%
%   Example
%     v1 = [2 0 0];
%     v2 = [0 3 0];
%     vectorCross3d(v1, v2)
%     ans =
%         0   0   6
%
%
%   Class support for inputs A,B:
%      float: double, single
%
%   See also DOT.

%   Sven Holcombe

% needed_colons = max([3, length(size(a)), length(size(b))]) - 3;
% tmp_colon = {':'};
% clnSet = tmp_colon(ones(1, needed_colons));
% 
% c = bsxfun(@times, a(:,[2 3 1],clnSet{:}), b(:,[3 1 2],clnSet{:})) - ...
%     bsxfun(@times, b(:,[2 3 1],clnSet{:}), a(:,[3 1 2],clnSet{:}));

sza = size(a);
szb = size(b);

% Initialise c to the size of a or b, whichever has more dimensions. If
% they have the same dimensions, initialise to the larger of the two
switch sign(numel(sza) - numel(szb))
    case 1
        c = zeros(sza);
    case -1
        c = zeros(szb);
    otherwise
        c = zeros(max(sza, szb));
end

c(:) =  bsxfun(@times, a(:,[2 3 1],:), b(:,[3 1 2],:)) - ...
        bsxfun(@times, b(:,[2 3 1],:), a(:,[3 1 2],:));
end
