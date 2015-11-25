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

function y = tracerPoints(position, directions, N, option) 
    %directions = tableau de vecteur
    %N fois tracerUneLigne
    for i = 1:N
        tracerUneLigne();
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
    positionCollision = collision(2);
    typeCollision = collision(1);
    
    %calculer distance
    distance = calculerdistance()
    
    if(typeCollision == 0)
        angle = calculerAngle2Vecteur(direction, normal);

        n = indiceRefraction(position, direction, option)
        critique = verifierAngleCritique(angle, n(1), n(2));
        if (critique == 0)
            nouvelleDirection = calculerDirectionRefraction(angle, n(1), n(2));
        else
            nouvelleDirection = calculerDirectionReflexion(angle, n(1), n(2));
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
    
end

function y = calculerAngle2Vecteur(Vecteur1, Vecteur2)
    CosTheta = dot(Vecteur1,Vecteur2)/(norm(Vecteur1)*norm(Vecteur2));
    y = acos(CosTheta);
end

function y = calculerDirectionRefraction()

end

function y = calculerDirectionReflexion()

end

function y = verifierAngleCritique(angle, n1, n2)
    %formule AngleCritique
    
    %1 = réflexion
    %0 = refraction
end