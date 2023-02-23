//MESHLESS METHOD OF FUNDAMENTAL SOLUTION (MFS)
//Autor: Guilherme Costa
//UFC - Laboratório de Hidráulica Computacional - jul.2021

//Description: Method of Fundamental Solutions basic equations and functions for Laplacian Operator in 2D

//Resets console
clear();
clc();
clf(0);
clf(1);

//#[1] MFS CORE

function phi = phiFS(x,xi,y,yi)
    phi = log(sqrt((x-xi)^2+(y-yi)^2))/(2*%pi)
endfunction

function line = phiLine(point,pointVec)
    for i=1:size(pointVec)(1)
        line(i) = phiFS(point(1),pointVec(i,1),point(2),pointVec(i,2))
    end
endfunction

function phi = phiDFS(x1,x2,y1,y2,n1,n2)
    phi = (-(x2-x1)/(2*%pi*((x2-x1)^2+(y2-y1)^2)))*n1+(-(y2-y1)/(2*%pi*((x2-x1)^2+(y2-y1)^2)))*n2
endfunction

function line = dphiLine(point,pointVec,n)
    for i=1:size(pointVec)(1)
        line(i) = phiDFS(point(1),pointVec(i,1),point(2),pointVec(i,2),n(1),n(2))
    end
endfunction

//#[2] MFS COMPUTATION

// Main Matrix
function mMatrix = mainMatrix(boundaryMatrix, offsetMatrix)
    for i=1:size(boundaryMatrix)(1)
        if boundaryMatrix(i,3)==1
            then
            mMatrix(i,:) = phiLine(boundaryMatrix(i,:),offsetMatrix)
        else
            if boundaryMatrix(i,3)==2
                then
                n(1) = boundaryMatrix(i,5)
                n(2) = boundaryMatrix(i,6)
                mMatrix(i,:) = dphiLine(boundaryMatrix(i,:),offsetMatrix,n(:))
            end
        end
    end
endfunction

function dMatrix = domainMainMatrix(domainMatrix, offsetMatrix)
    for i=1:size(domainMatrix)(1)
        dMatrix(i,:) = phiLine(domainMatrix(i,:), offsetMatrix)
    end
endfunction

function fgVector = boundaryConditionsVector(boundaryMatrix)
    for i=1:size(boundaryMatrix)(1)
        fgVector(i,1) = boundaryMatrix(i,4)
    end
endfunction

function numericalResult = fsmResult(offsetValue, domainMatrix, boundaryMatrix)
    offSetMatrix = getOffSetMatrix(offsetValue, boundaryMatrix)
    fgVector = boundaryConditionsVector(boundaryMatrix)
    
    mMatrix = mainMatrix(boundaryMatrix, offSetMatrix)
    invMMatrix = inv(mMatrix)
    lambda = invMMatrix * fgVector
    
    dMatrix = domainMainMatrix(domainMatrix, offSetMatrix)
    
    numericalResult = dMatrix * lambda
endfunction

