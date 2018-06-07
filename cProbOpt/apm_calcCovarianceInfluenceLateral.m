function [covInfluenceSys,covInfluenceRand] = apm_calcCovarianceInfluenceLateral(x,spots,ucm,expDose_ij,showWaitbar)

if nargin < 5
    showWaitbar = false;
end


nVox = numel(x);
nSpots = numel(spots);

covInfluenceSys = zeros(nVox,nSpots,nVox,nSpots);
covInfluenceRand = zeros(nVox,nSpots,nVox,nSpots);

for i = 1:nVox
    for l=i:nVox
        for j = 1:nSpots
            for m = 1:nSpots
                %Compute the element
                [covInfluenceSys(i,j,l,m) , covInfluenceRand(i,j,l,m)] = apm_calcCovInfluenceElement(i,j,l,m,x,spots,ucm,expDose_ij);
                
                %Symmetry
                covInfluenceSys(l,m,i,j) = covInfluenceSys(i,j,l,m);
                covInfluenceRand(l,m,i,j) = covInfluenceRand(i,j,l,m);
            end
        end
        
    end

    if (showWaitbar)
        if i == 1
            h = waitbar(i/nVox,'Computing covariance influence...');
        else
            waitbar(i/nVox,h,'Computing covariance influence...');
        end
    else
        matRad_progress(i,nVox);
    end
    
end

delete(h);

end

