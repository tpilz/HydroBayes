function [x , p_x ] = de_mc ( prior , pdf ,N ,T ,d)
% Differential Evolution Markov Chain (DE -MC) algorithm

gamma_RWM = 2.38/ sqrt (2* d ); % Calculate default jump rate
x = nan (T ,d ,N ); p_x = nan (T , N ); % Preallocate chains and density
X = prior (N ,d ); % Create initial population
for i = 1: N , p_X (i ,1) = pdf (X(i ,1: d )); end % Compute density initial population
    x (1 ,1: d ,1: N) = reshape (X',1 ,d ,N ); p_x (1 ,1: N) = p_X'; % Store initial states and density
    for i = 1: N , R (i ,1: N -1) = setdiff (1: N , i ); end % R- matrix : index of chains for DE
    for t = 2: T , % Dynamic part : Evolution of N chains 
        [~ , draw ] = sort ( rand (N -1 , N )); % Permute [1 ,... ,N -1] N times
        g = randsample ([ gamma_RWM 1] ,1 , true ,[0.9 0.1]); % Select gamma : 90/10 mix [ default 1] 
    for i = 1:N , % Create proposals and accept / reject
        a = R (i , draw (1 , i )); b = R(i , draw (2 , i )); % Extract a and b not equal i 
        Xp (i ,1: d) = X (i ,1: d) + g *( X(a ,1: d) -X(b ,1: d ))...
        + 1e-6* randn (1 , d ); % Create ith proposal with diff . evol . 
        p_Xp (i ,1) = pdf ( Xp (i ,1: d )); % Calculate density ith proposal
        p_acc = min (1 , p_Xp (i ,1)/ p_X (i ,1)); % Compute acceptance probability 
        if p_acc > rand , % p_acc larger than U [0 ,1]?
            X(i ,1: d) = Xp (i ,1: d ); p_X (i ,1) = p_Xp (i ,1); % True : Accept proposal 
    end
end 
x(t ,1: d ,1: N) = reshape (X',1 ,d ,N ); p_x (t ,1: N) = p_X'; % Append current X and density
end % End dynamic part 