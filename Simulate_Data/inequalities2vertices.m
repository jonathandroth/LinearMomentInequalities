function [V,nr] = inequalities2vertices(A,b)
% CON2VERT - convert a convex set of linear inequalities into the set
%            of vertices at the intersections of those inequalities;i.e.,
%            solve the "vertex enumeration" problem. Additionally,
%            identify redundant entries in the list of inequalities.
% 
% V = inequalities2vertices(A,b)
% [V,nr] = inequalities2vertices(A,b)
% 
% Converts the polyhedral set (a convex polygon, polyhedron, polytope that 
% may not be bounded) defined by the system of inequalities A*x <= b into 
% a list of vertices V. Each ROW of V is a vertex. For n variables:
% A = m x n matrix, where m >= n (m constraints, n variables)
%            (more precisely, rank(A)=n so that at least one vertex exists)
% b = m x 1 vector (m constraints)
% V = p x n matrix (p vertices, n variables)
% nr = list of the rows in A which are NOT redundant constraints
% 
% NOTES: (1) This program employs a primal-dual polytope method.
%        (2) In dimensions higher than 2, redundant vertices can
%            appear using this method. This program detects redundancies
%            at up to 6 digits of precision, then returns the
%            unique vertices.
%        (3) This program requires that the feasible region have some
%            finite extent in all dimensions. For example, the feasible
%            region cannot be a line segment in 2-D space, or a plane
%            in 3-D space.
%        (4) At least two dimensions are required.
%        (5) This program was written by Gregory Cox (National University 
%            of Singapore) in Feb 2020. It is based on the function 
%            con2vert written by Michael Klender in 2005. This function
%            generalizes con2vert by allowing for unbounded polyhedral sets
%            and handles the case m<n.  
%

c = A\b;
if ~all(A*c < b)
    fun = @(point) obj(point,A,b);
    [c,~,ef] = fminsearch(fun,c);
    if ef ~= 1
        error('Unable to locate a point within the interior of a feasible region.')
    end
end
b = b - A*c;
D = A ./ repmat(b,[1 size(A,2)]);
k = convhulln([D;zeros(1,size(D,2))]);
%[k,v1] = convhulln(D);
%if v2 > v1
%    error('Non-bounding constraints detected. (Consider box constraints on variables.)')
%end
k(any(k==size(D,1)+1,2),:) = []; % Remove all rows that contain the zeros point. This deals with unbounded set. 
nr = unique(k(:));
G  = zeros(size(k,1),size(D,2));
for ix = 1:size(k,1)
    F = D(k(ix,:),:);
    G(ix,:) = pinv(F)*ones(size(F,1),1);
end
G(any(~isfinite(G),2),:)=[]; % This deletes the NaN's and Inf's, which can only arise when inequalities are parallel and the set is not bounded in a direction. 
V = G + repmat(c',[size(G,1),1]);
[~,I]=unique(num2str(V,6),'rows');
V=V(I,:);
return
function dn = obj(cn,An,bn)
dn = An*cn-bn;
kn=(dn>=-1e-15);
dn(kn)=dn(kn)+1;
dn = max([0;dn]);
return

