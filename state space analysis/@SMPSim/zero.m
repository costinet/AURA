function [c,ceq] = zero(obj,t)
%ZERO This function creates the nonlinear function needed for fmincon to
%find a zero derivative within the time interval of a single state
%   I_solemnly_swear_that_I_am_up_to_no_good is a GLOBAL variable that is
%   created, altered, and cleared by VfwdIrev to iterate through the
%   different rows of the A and B matrix. This is because we are only
%   wanting to find when one state is equal to zero. The other stats can be
%   anything so they are not included in the constraints.

global  I_solemnly_swear_that_I_am_up_to_no_good;

i = I_solemnly_swear_that_I_am_up_to_no_good(1);
j = I_solemnly_swear_that_I_am_up_to_no_good(2);
k = I_solemnly_swear_that_I_am_up_to_no_good(3);

As = obj.As;
Bs = obj.Bs;
Xs = obj.Xs;
u = obj.u;

ceq_almost = As(:,:,k)*expm(As(:,:,k).*t)*Xs(:,j)+expm(As(:,:,k).*t)*Bs(:,:,k)*u;
ceq = ceq_almost(i);
c = [];

end
