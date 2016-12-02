%% Canonical Correlation Analysis on Riemannian Manifolds
%
% Canonical correlation analysis (CCA) is a widely used statistical technique 
% to capture correlations between two sets of multi-variate random variables.
% It has been applied to a multitude of applications in computer vision, 
% medical imaging and machine learning. 
% The classical formulation assumes that the data lives in a pair of vector 
% spaces which makes its use in certain important scientific domains problematic. 
% In this project, using the space of SPD matrices as a concrete example, 
% we give a principled generalization of the well known CCA to the 
% Riemannian setting. Our CCA algorithm operates on the product Riemannian 
% manifold representing SPD matrix-valued fields to identify meaningful 
% statistical relationships on the product Riemannian manifold. 
%
%
% * [ECCV2014] Hyunwoo J. Kim, Nagesh Adluru, Barbara B. Bendlin, Sterling C. Johnson, 
% Baba C. Vemuri, Vikas Singh, Canonical Correlation Analysis on Riemannian 
% Manifolds and its Applications, In European Conference on Computer Vision 
% (ECCV), September 2014.
%
% * [RCCV2015] Hyunwoo J. Kim, Nagesh Adluru, Barbara B. Bendlin, Sterling C. Johnson, 
% Baba C. Vemuri, Vikas Singh, Canonical Correlation Analysis on SPD(n) 
% Manifolds, Riemannian Computing and Statistical Inferences in Computer 
% Vision (RCCV), 2015.
%
% Project page:
% http://pages.cs.wisc.edu/~hwkim/projects/riem-cca/
%
% Github repository:
% http://github.com/MLman/riem-cca/ 
%
% Github page:
% http://mlman.github.io/riem-cca/
%
% The last update by <http://pages.cs.wisc.edu/~hwkim/index.html Hyunwoo J Kim>  2016/03/30 15:21:23 (CST)

%% Demos
%% Riemannian CCA and Euc CCA
clear;
close;
Ns =  [100 100 100 1000 1000 1000]; % number of samples
nexp = length(Ns);
exp_ids = 1:nexp;
results  = cell(nexp,1);
c1 = 1;
c2 = 5;
dimX = 3;
dimY = 3;




for iexp = 1:nexp
    % Generate two groups of samples in the space of X.
    N = Ns(iexp); % Number of samples
    [X, Y, group1, group2] = generate_data_for_demo1(N, dimX,dimY, c1,c2);   
    
    % CCA comparison
    % Riem-CCA-ga
    [A, B, r, U,V,stats] = rcca_ga(X, Y);
    
    % Euclidean CCA
    Y_euc = mxstack2mat(Y);
    X_euc = mxstack2mat(X);
    [Aeuc,Beuc,reuc,Ueuc,Veuc,stats_euc] = canoncorr(X_euc,Y_euc); 
    
    % Save resultss
    rcca.r = r(1);
    rcca.u = U(:,1);
    rcca.v = V(:,1);
    rcca.A = A(:,1);
    rcca.B = B(:,1);
    
    euc_cca.r = reuc(1);
    euc_cca.u = Ueuc(:,1);
    euc_cca.v = Veuc(:,1);
    euc_cca.A = Aeuc(:,1);
    euc_cca.B = Beuc(:,1);
    
    result = [];
    result.rcca = rcca;
    result.euc_cca = euc_cca;
    result.group1 = group1;
    result.group2 = group2;
    result.N = N;
    results{iexp} = result;
end 

%% Print experiment results (correlation)
fprintf('Exp ID: Correlation of Riem-CCA-ga, Euclidean CCA\n');
for iexp = 1:length(Ns)
    tmp =results{iexp};
    fprintf('Exp %d: Riem-CCA-ga %.4f, CCA %.4f\n',iexp, tmp.rcca.r, tmp.euc_cca.r);
end

%% Plot projections by Riem-CCA and Euclidean CCA
hFig = figure(1);
set(hFig, 'Position', [1 1 2400 1000])
for iexp = 1:length(Ns)
    % Plot
    tmp =results{iexp};
    
    subplot(2,nexp,iexp);
    plot(tmp.rcca.u,tmp.rcca.v,'o');
    hold on;
    title(['Riem-CCA-ga (corr=',sprintf('%.2f',tmp.rcca.r),')']);
    xlabel('Proj_w(X)');
    ylabel('Proj_w(Y)');
    hold off;
    % Plot
    subplot(2,nexp,nexp+iexp);
    plot(tmp.euc_cca.u,tmp.euc_cca.v,'o');
    hold on;
    title(['Euc-CCA (corr=',sprintf('%.2f',tmp.euc_cca.r),')']);
    xlabel('Proj_w(X)');
    ylabel('Proj_w(Y)');
    hold off;
end
