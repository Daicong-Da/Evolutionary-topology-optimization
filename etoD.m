%% EVOLUTIONARY TO WITH SMOOTH BOUNDARY REPRESENTATION
function etoD(nelx,nely,volfrac,er,rmin,ctp)
 %% INITIALIZATION
 vol = 1; change = 1;  ij = 0; xmin = 1e-6;
 vx = ones(nely,nelx);
 %% MATERIAL PROPERTIES
 E0 = 1; Emin = E0*xmin; nu = 0.3; pen = 3.0;
 %% PREPARE FINITE ELEMENT ANALYSIS
 A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
 A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
 B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
 B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
 KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
 nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+ nelx);
 edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
 edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
 iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
 jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
 %% ELEMENTAL NODES AND COORDINATES
 nodelast = reshape(nodenrs(1:end-1,1:end-1),nelx*nely,1);
 elenod = repmat(nodelast,1,4)+repmat([1 nely+[2 1] 0 ],nelx*nely,1);
 [nodex,nodey] = meshgrid(0:1:nelx,nely:-1:0);
 %% DEFINE LOADS AND SUPPORTS
switch(ctp)
 case 1 % HALF-MBB BEAM
 F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
 fixeddofs = union(1:2:2*(nely+1), 2*(nely+1)*(nelx+1));
 case 2 % CANTILEVER
 F = sparse(2*(nely+1)*(nelx+1)-nely,1,-1,2*(nely+1)*(nelx+1),1);
 fixeddofs = (1:2*(nely+1));
 case 3 % HALF-WHEEL
 F = sparse(2*(nely+1)*(nelx/2+1),1,-1,2*(nely+1)*(nelx+1),1);
 fixeddofs = union(2*nely+1:2*(nely+1), 2*(nely+1)*(nelx+1));
 end
 U = zeros(2*(nely+1)*(nelx+1),1);
 alldofs = (1:2*(nely+1)*(nelx+1));
 freedofs = setdiff(alldofs,fixeddofs);
 %% PREPARE FILTER
 iH = ones((nelx+1)*(nely+1)*(2*(ceil(rmin)+1))^2,1);
 jH = ones(size(iH)); sH = zeros(size(iH)); k = 0;
 [elex,eley] = meshgrid(1.5:nelx+0.5,1.5:nely+0.5);
 for i1 = 1:nelx+1
   for j1 = 1:nely+1
       e1 = (i1-1)*(nely+1)+j1;
         for i2 = max(i1-ceil(rmin),1):min(i1+ceil(rmin)-1,nelx)
             for j2 = max(j1-ceil(rmin),1):min(j1+ceil(rmin)-1,nely)
                 e2 = (i2-1)*nely+j2; k = k+1; iH(k) = e1;
                 jH(k) = e2;
                 sH(k) = max(0,rmin-sqrt((i1-elex(j2,i2))^2+(j1-eley(j2,i2))^2));
             end
         end
    end
 end
 H = sparse(iH,jH,sH); Hs = sum(H,2);
 %% START ITERATION
while change > 0.0001
 ij = ij + 1; vol = max(vol*(1-er), volfrac);
 if ij > 1; olddcnd = dcnd; end
 %% FE-ANALYSIS
 sK = reshape(KE(:)*(vx(:)'*E0+(1-vx(:))'*Emin),64*nelx*nely,1);
 K = sparse(iK,jK,sK); K = (K+K')/2;
 U(freedofs) = K(freedofs,freedofs)\F(freedofs);
 %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
 ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
 c(ij) = sum(sum((vx.*E0+(1-vx).*Emin).*ce));
 dc = ((1-vx)*xmin.^(pen-1)+vx)*E0.*ce;
 %% FILTERING/MODIFICATION OF NODAL SENSITIVITIES
 dcnd = reshape((H*dc(:)./Hs),nely+1,nelx+1);
 if ij > 1; dcnd = (dcnd+olddcnd)/2.; end
 %% OPTIMALITY CRITERIA UPDATE OF ELEMENT VOLUME
 l1 = min(dcnd(:)); l2 = max(dcnd(:));
 while (l2-l1)/abs(l1+l2) > 1.0e-9
 ls = (l1+l2)/2.0;
 dcth = dcnd(:)-ls;
   for i = 1:nely*nelx
       if min(dcth(elenod(i,:))) > 0
           vx(i) = 1;
       elseif max(dcth(elenod(i,:))) < 0
           vx(i) = 0; 
       else
           ngrid = 40; [s,t] = meshgrid(-1+1/ngrid:2/ngrid:1-1/ngrid,-1+1/ngrid:2/ngrid:1-1/ngrid);
           ps = (1 - s(:)).*(1 - t(:))/4 * dcth(elenod(i,1)) + (1 + s(:)).*(1 - t(:))/4 * dcth(elenod(i,2))...
               + (1 + s(:)).*(1 + t(:))/4 * dcth(elenod(i,3)) + (1-s(:)).*(1 + t(:))/4 * dcth(elenod(i,4));
           vx(i) = length(find( ps >= 0 ))/length(s(:));
       end
   end
 if mean(vx(:)) - vol > 0;
 l1 = ls;
 else
 l2 = ls;
 end
 end
 %% PLOT RESULTS
 if ij > 10; change =abs(sum(c(ij-9:ij-5))-sum(c(ij-4:ij)))/sum(c(ij-4:ij)); end
 fprintf('It.:%3i Obj.:%8.4f Vol.:%4.3f ch.:%4.5f\n',(ij),c(ij),mean(vx(:)),change);
 contourf(nodex, nodey, (dcnd-ls), [0 0] ); axis equal;  axis tight; axis off; pause(1e-6);
end