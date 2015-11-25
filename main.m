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

function y = positionBloc()
    
end

function y = indiceRefraction(position, option)
    
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

        n1 = indiceRefraction(position, option)
        critique = verifierAngleCritique(angle, n1, n2);
        if (critique == 0)
            nouvelleDirection = calculerDirectionRefraction(angle, n1, n2);
        else
            nouvelleDirection = calculerDirectionReflexion(angle, n1, n2);
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