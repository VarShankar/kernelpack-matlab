classdef OpProperties
    %OPPROPERTIES KernelPack-style operator metadata.

    properties
        selectdim (1,1) double = 0
        decompose (1,1) logical = true
        storeWeights (1,1) logical = true
        recordStencils (1,1) logical = false
        nosolve (1,1) logical = false
        OverlapLoad (1,1) double = 0.5
        UseParallel (1,1) logical = false
    end

    methods
        function obj = OpProperties(varargin)
            if nargin == 0
                return;
            end
            parser = inputParser();
            parser.addParameter('selectdim', obj.selectdim);
            parser.addParameter('decompose', obj.decompose);
            parser.addParameter('storeWeights', obj.storeWeights);
            parser.addParameter('recordStencils', obj.recordStencils);
            parser.addParameter('nosolve', obj.nosolve);
            parser.addParameter('OverlapLoad', obj.OverlapLoad);
            parser.addParameter('UseParallel', obj.UseParallel);
            parser.parse(varargin{:});

            fields = fieldnames(parser.Results);
            for k = 1:numel(fields)
                obj.(fields{k}) = parser.Results.(fields{k});
            end
        end
    end
end
