function main
    %little matlab fairy with great magic
    nbreDePoints = 1024;
    positionInit = [-10 -10 15];
    % generer les directions
    largeur = sqrt(7^2 + 7^2) + 1;
    hauteur = 20 + 1;
    coinBloc = [-0.5 -0.5 5];    
    directions = zeros(3);
    count = 1;
    for i = 1:sqrt(nbreDePoints)
        for j = 1:sqrt(nbreDePoints)
            directions(count, :) =  [coinBloc(1)+largeur/i coinBloc(2)+hauteur/j coinBloc(3)] - positionInit;
            count = count + 1;
        end
    end
    option = 1;
    tracerPoints(positionInit, directions, nbreDePoints, option);
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
        0 -1 0;
        1 0 0;
        0 1 0;
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
        hold on
        detailPoint = tracerUneLigne(position, directions(i), option);
        if(detailPoint(1) > 1)
            point = detailPoint(2)*directions(i) + position;
            scatter3(point(1),point(2),point(3)); %TODO : add color LOL
        end
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
        
        pointCollision = tracerUneLigne(positionCollision, nouvelleDirection, option);
    else
        pointCollision = [typeCollision distance];
    end
    
    y = [pointCollision(1) (pointCollision(2)+distance)];
end

function y = calculerdistance(position1, position2)
    X = [position1; position2] ;
    y = pdist(X, 'euclidean');
end

function y = calculerCollision(position, direction)
    
    Bloc = positionBloc();
    petitBloc = positionPetitBloc();
    grosseBoite = positionGrosseBoite();
    vraiIntersection = [0 0 0];
    minDistance = 1000000000;
    normalSurface = [0 0 0];
    typeCollision = 0;
    % check des collisions avec les face laterales des deux blocs
    for i = 1:4
        intersection = intersectLinePlane([position position+direction], [Bloc(i,:) Bloc(i+1, :) Bloc(i+4, :)]);

        %intersection avec le bloc de verre (ou de polymere bizzare...)
        if(isBetweenTwoPoints(intersection, Bloc(mod(i+1,4), :), Bloc(i+4,:)) == 1)
            temp = intersection - position;
            j=1;
            while(temp(j) ==0 && j < 4)
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

        intersectionPetit = intersectLinePlane([position position+direction], [petitBloc(i,:) petitBloc(i+1, :) petitBloc(i+4, :)]);
        %intersection avec le bloc d'acier
        if(isBetweenTwoPoints(intersection, petitBloc(mod(i+1,4), :), petitBloc(i+4,:)) == 1)
            temp = intersectionPetit - position;
            j=1;
            while(temp(j) ==0 && j < 4)
                j = j+1;
            end

            a = temp(j)/direction(j);
            if(a > 0 && norm(temp) < minDistance)
                minDistance = norm(temp);
                vraiIntersection = intersection;
                normalTemp = normaleSurface();
                normalSurface = normalTemp(i, :);
                typeCollision = i + 1;

                %x
                % garder la surface la plus proche (plus petit x positif)
            end
        end
    end

    % intersection pour les faces superieurs
    intersection = intersectLinePlane([position position+direction], [Bloc(1,:) Bloc(2, :) Bloc(3, :)]);

    %intersection avec le bloc de verre (ou de polymere bizzare...)
    if(isBetweenTwoPoints(intersection, Bloc(1, :), Bloc(3,:)) == 1)
        temp = intersection - position;
        j=1;
        while(temp(j) ==0 && j < 4)
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

    intersectionPetit = intersectLinePlane([position position+direction], [petitBloc(1,:) petitBloc(2, :) petitBloc(3, :)]);
    %intersection avec le bloc d'acier
    if(isBetweenTwoPoints(intersection, petitBloc(1, :), petitBloc(3,:)) == 1)
        temp = intersectionPetit - position;
        j=1;
        while(temp(j) ==0 && j < 4)
            j = j+1;
        end

        a = temp(j)/direction(j);
        if(a > 0 && norm(temp) < minDistance)
            minDistance = norm(temp);
            vraiIntersection = intersection;
            normalTemp = normaleSurface();
            normalSurface = normalTemp(5, :);
            typeCollision = 6;
            %x
            % garder la surface la plus proche (plus petit x positif)
        end
    end


    % intersection pour les faces inferieurs
    intersection = intersectLinePlane([position position+direction], [Bloc(5,:) Bloc(6, :) Bloc(7, :)]);

    %intersection avec le bloc de verre (ou de polymere bizzare...)
    if(isBetweenTwoPoints(intersection, Bloc(5, :), Bloc(7,:)) == 1)
        temp = intersection - position;
        j=1;
        while(temp(j) ==0 && j < 4)
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

    intersectionPetit = intersectLinePlane([position position+direction], [petitBloc(5,:) petitBloc(6, :) petitBloc(7, :)]);
    %intersection avec le bloc d'acier
    if(isBetweenTwoPoints(intersection, petitBloc(1, :), petitBloc(3,:)) == 1)
        temp = intersectionPetit - position;
        j=1;
        while(temp(j) ==0 && j < 4)
            j = j+1;
        end

        a = temp(j)/direction(j);
        if(a > 0 && norm(temp) < minDistance)
            minDistance = norm(temp);
            vraiIntersection = intersection;
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