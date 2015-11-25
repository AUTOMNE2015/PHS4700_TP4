function main
    %little matlab fairy with great magic
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
    y = [1 0 0;
        -1 0 0;
        0 1 0;
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
             
    point1 = posBloc(1, :);
    point2 = posBloc(7, :);
    %n1
    if(option == 1)
        n1 = 1.0;
    else
        n1 = 1.33;
    end
    if(point1(1) < posArriere(1) && point2(1) > posArriere(1))
        if(point1(2) < posArriere(2) && point2(2) > posArriere(2))
             if(point1(3) < posArriere(3) && point2(3) > posArriere(3))
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
    if(point1(1) < posDevant(1) && point2(1) > posDevant(1))
        if(point1(2) < posDevant(2) && point2(2) > posDevant(2))
             if(point1(3) < posDevant(3) && point2(3) > posDevant(3))
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
    if(point1(1) <= point3(1) && point2(1) >= point3(1))
        if(point1(2) <= point3(2) && point2(2) >= point3(2))
             if(point1(3) <= point3(3) && point2(3) >= point3(3))
                 y = 1;
             end
        end
    end          
end


function y = tracerPoints(position, directions, nombrePoints, option) 
    %directions = tableau de vecteur
    %N fois tracerUneLigne
    for i = 1:nombrePoints
        tracerUneLigne(position, directions(i), option);
    end
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
    positionCollision = collision(2);
    normal = collision(3);
    
    %calculer distance
    distance = calculerdistance(position, positionCollision);
    
    if(typeCollision == 0)
        angle = calculerAngle2Vecteur(direction, normal);

        n = indiceRefraction(position, direction, option);
        critique = verifierAngleCritique(angle, n(1), n(2));
        %critique = 1 = réflexion
        %critique = 0 = refraction
        if (critique == 0)
            nouvelleDirection = calculerDirectionRefraction(direction, normal, n(1), n(2));
        else
            nouvelleDirection = calculerDirectionReflexion(direction, normal);
        end
        
        r = tracerUneLigne(positionCollision, nouvelleDirection, option);
    else
        r = [typeCollision distance];
    end
    
    y = [r(1) (r(2)+distance)];
end

function y = calculerdistance(position1, position2)
    X = [position1; position2] ;
    y = pdist(X, 'euclidean');
end

function y = calculerCollision(position, direction)
    
    Bloc = positionBloc();
    petitBloc = positionPetitBloc();
    grosseBoite = positionGrosseBoite();
    for i = 1:4
        intersection = intersectLinePlane([position position+direction], [Bloc(i,:) Bloc(i+1, :) Bloc(i+4, :)]);
        if(isBetweenTwoPoints(intersection, Bloc(mod(i+1,4), :), Bloc(i+4,:)) == 1)
            temp = position1 - position2;
            x = temp(1)/direction(1);
            if(x*direction == temp)
                %x
                % garder la surface la plus proche (plus petit x positif)
            end
        end
    end
    % refaire la verification avec 2 faces restantes
    %trouver le plus petit x > 0
    % retourner la surface
    typeCollision = 0;
    positionCollision = [0 0 0];
    normalSurface = [0 1 0];
    y = [typeCollision, positionCollision, normalSurface];
end

function y = calculerAngle2Vecteur(Vecteur1, Vecteur2)
    CosTheta = dot(Vecteur1,Vecteur2)/(norm(Vecteur1)*norm(Vecteur2));
    y = acos(CosTheta);
end

function y = calculerDirectionRefraction(direction, normale, n1, n2)
    % voir slide p48
    j = cross(direction, normale)/norm(cross(direction, normale));
    k = cross(normale, j); 
    s = (n1/n2)*(dot(direction,k)); 
    dir = -normale * sqrt(1-s^2) + k*s;
    y = dir/norm(dir); %unitaire
end

function y = calculerDirectionReflexion(direction, normale)
    %voir slide p.35
    res = directio - 2*dot(direction, normale)*normale; 
    y = res/norm(res);
end

function y = verifierAngleCritique(angle, n1, n2)
    %formule AngleCritique
    sb = sin(angle) * n1 / n2;
    if(sb > 1)
        y = 1;
    else
        y = 0;  
    end
    %1 = réflexion
    %0 = refraction
end

function point = intersectLinePlane(line, plane)

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

    % set indices of line and plane which are parallel to NaN
    par = abs(denom) < tol;
    point(par,:) = NaN;
end