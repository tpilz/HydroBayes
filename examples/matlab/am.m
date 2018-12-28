function [x , p_x ] = am ( prior , pdf ,T , d)
% Adaptive Metropolis (AM) algorithm

q = @ (C ,d ) mvnrnd ( zeros (1 , d), C ); % d- variate normal proposal distribution
C = (2.38/ sqrt (d ))^2 * eye (d ); % Covariance matrix proposal distribution 
x = nan (T , d ); p_x = nan (T ,1); % Preallocate memory for chain and density
x (1 ,1: d ) = prior (1 , d ); % Initialize chain by sampling from prior 
p_x (1) = pdf (x (1 ,1: d )); % Compute density initial state chain

for t = 2: T , % Dynamic part : evolution of chain
% --------------- Adaptation covariance matrix of proposal distribution ---------------- 
if ( mod (t ,10) == 0 )
C = (2.38/ sqrt (d ))^2 * ( cov (x (1: t -1 ,1: d )) + 1e-4* eye (d )); 
% Note : recursive formulae for C much more CPU - efficient !
end 
% --------------------------------- End adaptation -------------------------------------
xp = x(t -1 ,1: d) + q (C ,d ); % Create candidate point 
p_xp = pdf ( xp ); % Calculate density proposal
p_acc = min (1 , p_xp / p_x (t -1)); % Compute p_accept 
if p_acc > rand , % p_acc larger than U [0 ,1]?
x(t ,1: d ) = xp ; p_x (t) = p_xp ; % True : accept proposal 
else
x(t ,1: d ) = x(t -1 ,1: d ); p_x ( t) = p_x (t -1); % False : copy old state 
end
end % End dynamic part