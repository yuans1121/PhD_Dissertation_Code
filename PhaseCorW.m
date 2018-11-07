function PhaseCorrW=PhaseCorW(y,PhaseTau,lambda,theta)
% Correction term for the phase

global My Mx MaskX

Y=ones(My,Mx);
for ii=1:1:Mx
    Y(:,ii)=y; 
end
Y=flipud(Y);
    
PhaseCorr=PhaseTau-(2*pi/lambda)*Y*sin(theta);
PhaseCorrW=wrapToPi(PhaseCorr);

end