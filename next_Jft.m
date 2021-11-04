
%This function takes a previous period's offerings Jftminus1, and the
%transition probs p_o and p_d. 
%It generates the next periods J_ft using the given transition
%probabilities

function Jft = next_Jft( Jftminus1 , p_o, p_d )

draws = rand( size(Jftminus1) );

%If a product is out in Jftminus1, it is in with probability p_o
%If a product is in in Jftminus1, it is in with probability 1 - p_d
Jft = (Jftminus1 == 0) .* (draws < p_o) + (Jftminus1 == 1) .* (draws < (1-p_d) );


end
