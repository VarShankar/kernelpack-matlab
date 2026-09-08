function study = moving_surface_adr_tp_literature_comparison(varargin)
%MOVING_SURFACE_ADR_TP_LITERATURE_COMPARISON Final literature comparison.
%   Compatibility wrapper for the systematic xi=6 threshold sweep.  The
%   comparison intentionally asks for the coarsest final-method point cloud
%   that matches each published manufactured-solution error threshold.

study = moving_surface_adr_tp_literature_threshold_sweep(varargin{:});
end
