function [x , p_x ] = par_dream ( prior , pdf ,N ,T , d)
% DiffeRential Evolution Adaptive Metropolis ( DREAM ) algorithm

[ delta ,c , c_star , n_CR , p_g ] = deal (3 ,0.1 ,1e-12 ,3 ,0.2); % Default of algorithmic parameters
x = nan (T ,d ,N ); p_x = nan (T , N ); % Preallocate chains and density 
CR = [1: n_CR ]/ n_CR ; p_CR = ones (1 , n_CR )/ n_CR ; % Crossover values and select . prob .
[J , n_id ] = deal ( zeros (1 , n_CR )); % Variables selection prob . crossover 
for i = 1: N , R (i ,1: N -1) = setdiff (1: N , i ); end % R- matrix : index of chains for DE

X = prior (N ,d ); % Create initial population
for i = 1: N , p_X (i ,1) = pdf (X(i ,1: d )); end % Compute density initial population 
x (1 ,1: d ,1: N) = reshape (X',1 ,d ,N ); p_x (1 ,1: N) = p_X'; % Store initial states and density

for t = 2: T , % Dynamic part : Evolution of N chains
    [~ , draw ] = sort ( rand (N -1 , N )); % Permute [1 ,... ,N -1] N times 
    dX = zeros (N , d ); % Set N jump vectors to zero
    id = zeros(1,N);
    lambda = unifrnd (-c ,c ,N ,1); % Draw N lambda values 
    std_X = std (X ); % Compute std each sampling dimension
    for i = 1:N , % Create proposals 
        D = randsample ([1: delta ] ,1 , 'true'); % Select delta ( equal select . prob .)
        a = R (i , draw (1: D ,i )); b = R(i , draw (D +1:2* D , i )); % Extract vectors a and b not equal i 
        id (i) = randsample (1: n_CR ,1 , 'true' , p_CR ); % Select index of crossover value
        idt = id(i);
        z = rand (1 , d ); % Draw d values from U[0 ,1] 
        A = find (z < CR (idt)); % Derive subset A selected dimensions
        d_star = numel (A ); % How many dimensions sampled ? 
        if ( d_star == 0) , [~ , A] = min (z ); d_star = 1; end % A must contain at least one value
        gamma_d = 2.38/ sqrt (2* D* d_star ); % Calculate jump rate 
        g = randsample ([ gamma_d 1] ,1 , 'true' ,[1 - p_g p_g ]); % Select gamma : 80/20 mix [ default 1]
        dX (i ,A ) = c_star * randn (1 , d_star ) + ... 
        (1+ lambda ( i ))* g * sum (X(a ,A) -X(b , A ) ,1); % Compute ith jump diff . evol .
        Xp (i ,1: d) = X (i ,1: d) + dX (i ,1: d ); % Compute ith proposal 
    end
    parfor i = 1:N , % Accept / reject proposals ( parallel ) 
        p_Xp (i ,1) = pdf ( Xp (i ,1: d )); % Calculate density ith proposal
        p_acc = min (1 , p_Xp (i ,1)./ p_X (i ,1)); % Compute acceptance probability 
        if p_acc > rand , % p_acc larger than U [0 ,1]?
            X(i ,:) = Xp (i ,1: d ); p_X (i ,1) = p_Xp (i ,1); % True : Accept proposal 
        else
            dX (i,:) = zeros(1,d); % Set jump back to zero for p_CR 
        end     
    end 
    for i = i:N,
        idt = id(i);
        J( idt) = J(idt ) + sum (( dX (i ,1: d )./ std_X ).^2); % Update jump distance crossover idx 
        n_id(idt) = n_id(idt) + 1; % How many times idx crossover used
    end
    x(t ,1: d ,1: N) = reshape (X',1 ,d ,N ); p_x (t ,1: N) = p_X'; % Append current X and density
    %if t < T /10 , p_CR = J ./ n_id ; p_CR = p_CR / sum ( p_CR ); end % Update selection prob . crossover 
    %[X , p_X ] = check (X , mean ( log ( p_x ( ceil (t /2): t ,1: N )))); % Check for outlier and correct
end % End dynamic part 