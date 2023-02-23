//METHOD OF FUNDAMENTAL SOLUTIONS (MFS)
//Autor: Guilherme Costa
//Contato: guilherme3218@yahoo.com.br
//UFC - Laboratório de Hidráulica Computacional - ago.2019

//Descrição: Rotina para gerar dados da aplicação do FSM sob diferentes critérios de variação de parâmetros
//Utilizando funcoes analiticas para equacao de Laplace

//Reseta o console
clear();
clc();
clf(0);
clf(1);

//#[1] Importacao de pontos do contorno

/*Modificar o caminho do arquivo para o arquivo correspondente
Observar formatacao:
1- O arquivo será do tipo xls, arquivo do excel
2- Cada coluna correspondera a um tipo de informacao:
    [1] coordenadas x
    [2] coordenadas y
    [3] tipo de condicao de contorno
    [4] valor da condicao de contorno
    [5] coordenada x do vetor unitario normal ao ponto
    [6] coordenada y do vetor unitario normal ao ponto
*/

sheet = readxls("C:\...\Database.xls")
s1 = sheet(1)
boundaryMatrix(:,1) = s1(:,1) //coordenadas x
boundaryMatrix(:,2) = s1(:,2) //coordenadas y
boundaryMatrix(:,3) = s1(:,3) //tipo de condicao de contorno (1 ou 2)
boundaryMatrix(:,4) = s1(:,4) //valor da condicao de contorno
boundaryMatrix(:,5) = s1(:,5) //coordenada x do vetor unitario normal ao ponto de contorno
boundaryMatrix(:,6) = s1(:,6) //coordenada y do vetor unitario normal ao ponto de contorno

//#[2] Funcoes para a criacao das matrizes principais do metodo

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

//#[3] Funcoes para obtencao dos vetores normais aos pontos do contorno

function normalAnsVector = getNormalVectors(boundaryMatrix)
    n = size(boundaryMatrix)(1)
    //#[3.1] Calculo dos vetores unitarios
    for i=1:n
        /*
        temp[1] - modulo dos vetores
        temp[2] - coordenada x trasportada para origem
        temp[3] - coordenada y trasnportada para origem
        temp[4] - coordenada x do vetor unitario
        temp[5] - coordenada y do vetor unitario
        */
        if i~=n
            then
            temp(i,1) = (((boundaryMatrix(i,1)-boundaryMatrix(i+1,1))^2)+((boundaryMatrix(i,2)-boundaryMatrix(i+1,2))^2))^(1/2)
            temp(i,2) = boundaryMatrix(i+1,1) - boundaryMatrix(i,1)
            temp(i,3) = boundaryMatrix(i+1,2) - boundaryMatrix(i,2)
        else
            temp(i,1) = (((boundaryMatrix(i,1)-boundaryMatrix(1,1))^2)+((boundaryMatrix(i,2)-boundaryMatrix(1,2))^2))^(1/2)
            temp(i,2) = boundaryMatrix(1,1) - boundaryMatrix(i,1)
            temp(i,3) = boundaryMatrix(1,2) - boundaryMatrix(i,2)
        end
        temp(i,4) = temp(i,2)/temp(i,1)
        temp(i,5) = temp(i,3)/temp(i,1)
    end
    //#[3.2] Determinacao dos vetores normais aos pontos, com base em seus vizinhos
    //anterior e posterior
    for i=1:n
        if i~=1
            then
            normalVector(i,1) = temp(i-1,4) - temp(i,4)
            normalVector(i,2) = temp(i-1,5) - temp(i,5)
        else
            normalVector(i,1) = temp(n,4) - temp(i,4)
            normalVector(i,2) = temp(n,5) - temp(i,5)
        end
        
        //Pontos colineares
        if normalVector(i,1)==0 & normalVector(i,2)==0
            then
            angle = acos(temp(i,4))
            normalVector(i,1) = cos(angle+%pi/2)
            normalVector(i,2) = sin(angle+%pi/2)
        end
        
        module = ((normalVector(i,1)^2)+(normalVector(i,2)^2))^(1/2)
        normalAnsVector(i,1) = normalVector(i,1)/module
        normalAnsVector(i,2) = normalVector(i,2)/module
    end
       //Pontos colineares e alinhados a (0,-1) (tratamento de exceção)
    for i=1:n
        if i~=1
            then
            if(temp(i,5)==-1 & temp(i-1,5)==-1)
                then
                normalAnsVector(i,1) = cos(2*%pi)
                normalAnsVector(i,2) = 0
            end
        else
            if(temp(i,5)==-1 & temp(n,5)==-1)
                then
                normalAnsVector(i,1) = cos(2*%pi)
                normalAnsVector(i,2) = 0
            end
        end
    end
endfunction

function centroidCoordinates = getCentroidCoordinates(boundaryMatrix)
    n = size(boundaryMatrix)(1)
    tempCoordinates(1) = 0
    tempCoordinates(2) = 0
    
    for i=1:n
        tempCoordinates(1) = tempCoordinates(1) + boundaryMatrix(i,1)
        tempCoordinates(2) = tempCoordinates(2) + boundaryMatrix(i,2)
    end
    
    centroidCoordinates(1,1) = tempCoordinates(1)/n
    centroidCoordinates(1,2) = tempCoordinates(2)/n
endfunction

function offMatrix = getOffSetMatrix(offset, boundaryMatrix)
    n = size(boundaryMatrix)(1)
    centroidCoordinates = getCentroidCoordinates(boundaryMatrix)
    normalVector = getNormalVectors(boundaryMatrix)
    
    for i=1:n
        offMatrix(i,1) = boundaryMatrix(i,1) + offset*normalVector(i,1)
        offMatrix(i,2) = boundaryMatrix(i,2) + offset*normalVector(i,2)
        
        distB = (((boundaryMatrix(i,1)-centroidCoordinates(1,1))^2)+((boundaryMatrix(i,2)-centroidCoordinates(1,2))^2))
        distO = (((offMatrix(i,1)-centroidCoordinates(1,1))^2)+((offMatrix(i,2)-centroidCoordinates(1,2))^2))
        
        if distB > distO
            then
            offMatrix(i,1) = boundaryMatrix(i,1) - offset*normalVector(i,1)
            offMatrix(i,2) = boundaryMatrix(i,2) - offset*normalVector(i,2)
        end
    end
endfunction

//#[4] Funcoes para determinacao dos pontos do dominio

function lengthValue = lengthFunction(boundaryMatrix)
    xMax = 0
    xMin = 0
    yMax = 0
    yMin = 0
    
    for i=1:size(boundaryMatrix)(1)
        if boundaryMatrix(i,1) > xMax
            xMax = boundaryMatrix(i,1)
            xMin = xMax
        end
    end
    
    for i=1:size(boundaryMatrix)(1)
        if boundaryMatrix(i,1) < xMin
            xMin = boundaryMatrix(i,1)
        end
    end
    
    for i=1:size(boundaryMatrix)(1)
        if boundaryMatrix(i,2) > yMax
            yMax = boundaryMatrix(i,2)
            yMin = yMax
        end
    end
    
    for i=1:size(boundaryMatrix)(1)
        if boundaryMatrix(i,2) < yMin
            yMin = boundaryMatrix(i,2)
        end
    end
    
    lengthValue(1) = xMax
    lengthValue(2) = xMin
    lengthValue(3) = yMax
    lengthValue(4) = yMin
    
endfunction

function orientation = whichOrientation(p,q,r)
    value = (q(2)-p(2))*(r(1)-q(1))-(r(2)-q(2))*(q(1)-p(1))
    
    if value == 0 
        then
        orientation = 0
    else
        if value > 0
            orientation = 1
        else
            if value < 0
                orientation = 2
            end
        end
    end
endfunction

function intersection = doIntersect(seg1, seg2)
    o1 = whichOrientation(seg1(1,:),seg1(2,:),seg2(1,:))
    o2 = whichOrientation(seg1(1,:),seg1(2,:),seg2(2,:))
    o3 = whichOrientation(seg2(1,:),seg2(2,:),seg1(1,:))
    o4 = whichOrientation(seg2(1,:),seg2(2,:),seg1(2,:))
    
    if (o1~=o2 & o3~=o4) 
        then
        intersection = 1
    else
        if (o1==0 | o2==0 | o3==0 | o4==0)
            intersection = 2
        else
            intersection = 2
        end
    end
endfunction

function yesOrNo = pointIsInside(point, boundarySegments, xMaxExtended)
    seg1(1,1) = point(1)
    seg1(1,2) = point(2)
    seg1(2,1) = xMaxExtended
    seg1(2,2) = point(2)
    
    count = 0
    for i=1:size(boundarySegments)(1)
        seg2(1,1) = boundarySegments(i,1)
        seg2(1,2) = boundarySegments(i,2)
        seg2(2,1) = boundarySegments(i,3)
        seg2(2,2) = boundarySegments(i,4)
        
        if doIntersect(seg1(:,:), seg2(:,:))==1
            then
            count = count+1
        end
    end

    if ((count/2)-round(count/2))==0 then
        yesOrNo = 2
    else
        if count == 0
            then
            yesOrNo = 2
        else
            yesOrNo = 1
        end
    end
endfunction

function domain = generateInsideDomain(deltaDomain, boundaryMatrix)
    lengthValue = lengthFunction(boundaryMatrix)
    
    xStep = linspace(lengthValue(1),lengthValue(2),deltaDomain)
    yStep = linspace(lengthValue(3),lengthValue(4),deltaDomain)
    
    count1 = 1
    for i=1:deltaDomain
        for j=1:deltaDomain
            generalDomain(count1,1) = xStep(i)
            generalDomain(count1,2) = yStep(j)
            count1 = count1+1
        end
    end
    
    for i=1:size(boundaryMatrix)(1)
        if i~=size(boundaryMatrix)(1)
            then
            bSegments(i,1) = boundaryMatrix(i,1)
            bSegments(i,2) = boundaryMatrix(i,2)
            bSegments(i,3) = boundaryMatrix(i+1,1)
            bSegments(i,4) = boundaryMatrix(i+1,2)
        else
            bSegments(i,1) = boundaryMatrix(i,1)
            bSegments(i,2) = boundaryMatrix(i,2)
            bSegments(i,3) = boundaryMatrix(1,1)
            bSegments(i,4) = boundaryMatrix(1,2)
        end
    end
    
    count2 = 1
    for i=1:size(generalDomain)(1)
        xMaxExtended = lengthValue(2)+1.5*(lengthValue(1)-lengthValue(2))
        yOn = pointIsInside(generalDomain(i,:), bSegments(:,:), xMaxExtended)
        if yOn == 1
            then
            domain(count2,1) = generalDomain(i,1)
            domain(count2,2) = generalDomain(i,2)
            count2 = count2+1
        end
    end
endfunction

//#[5] Funcoes da aplicacao do metodo

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

function numericalResult = fsmResult(offsetValue, deltaDomain, boundaryMatrix)
    offSetMatrix = getOffSetMatrix(offsetValue, boundaryMatrix)
    domainMatrix = generateInsideDomain(deltaDomain, boundaryMatrix)
    fgVector = boundaryConditionsVector(boundaryMatrix)
    
    mMatrix = mainMatrix(boundaryMatrix, offSetMatrix)
    invMMatrix = inv(mMatrix)
    lambda = invMMatrix * fgVector
    
    dMatrix = domainMainMatrix(boundaryMatrix, offSetMatrix)
    
    numericalResult = dMatrix * lambda
endfunction

//#[6] Funções de analise de resultados

function analyticalResult = analyticResult(deltaDomain, boundaryMatrix)
    domainMatrix = generateInsideDomain(deltaDomain, boundaryMatrix)
    
    for i=1:size(domainMatrix)(1)
        x = domainMatrix(i,1)
        y = domainMatrix(i,2)
        analyticalResult(i) = 0//aplica soluçao analitica 
    end
endfunction

function erroAns = erro(vector1, vector2)
    n = size(vector1)
    m = size(vector2)
    
    if n~=m
        then
        disp("Erro: the two vectors must have the same size")
        abort
    end
    
    for i=1:n
        erroAns(i) = vector1(i) - vector2(i)
    end
endfunction

function rmseAns = rmse(erroVector)
    temp = 0
    for i=1:size(erroVector)
        temp = temp + erroVector(i)^2
    end
    rmseAns = (temp/size(erroVector))^(1/2)
endfunction

//#[FINAL] TESTE
offsetValue = 1
deltaDomain = 1
scatter(boundaryMatrix(:,1),boundaryMatrix(:,2), 50, "filled diamond")
offsetPoints = getOffSetMatrix(offsetValue, boundaryMatrix)
scatter(offsetPoints(:,1),offsetPoints(:,2), "x")
insideDomain = generateInsideDomain(deltaDomain,boundaryMatrix)
scatter(insideDomain(:,1),insideDomain(:,2),1,".")
result = fsmResult(offsetValue, deltaDomain, boundaryMatrix)
