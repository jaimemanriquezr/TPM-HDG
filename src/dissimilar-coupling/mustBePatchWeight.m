function mustBePatchWeight(weight)
% Validates the patch_weight option of the coupled solvers:
% "jacobian" (exact ruled-surface area element) or "rectangular"
% (strip weight, valid for near-normal transfer paths).
mustBeMember(weight, ["jacobian", "rectangular"])
end
