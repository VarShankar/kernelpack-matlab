function out = pusl_advection_diffusion_reaction_common()
%PUSL_ADVECTION_DIFFUSION_REACTION_COMMON Shared manufactured problem for ADR PU-SL solvers.

base = pusl_advection_diffusion_common();

out = base;
out.reaction_rate = 0.35;
out.reaction = @reaction;
out.reaction_exact = @reaction_exact;
out.forcing = @forcing_with_reaction;
end

function values = reaction(t, state, X) %#ok<INUSD>
common = pusl_advection_diffusion_common();
rate = 0.35;
values = rate * (state(:) - 1.0);
expected = reaction_exact(t, X);
if numel(values) ~= numel(expected)
    error('kp:examples:BadReactionSize', 'Reaction callback returned the wrong size.');
end
end

function values = reaction_exact(t, X)
common = pusl_advection_diffusion_common();
rate = 0.35;
values = rate * exp(-t) .* common.phi(X);
end

function values = forcing_with_reaction(domain, nu, t, X)
common = pusl_advection_diffusion_common();
values = common.forcing(domain, nu, t, X) - reaction_exact(t, X);
end
