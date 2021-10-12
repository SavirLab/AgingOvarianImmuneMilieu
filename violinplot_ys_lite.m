function violins = violinplot_ys_lite(data, x,median_flag, varargin)
%Violinplots plots violin plots of some data and categories
%   VIOLINPLOT(DATA) plots a violin of a double vector DATA
%
%   VIOLINPLOT(DATAMATRIX) plots violins for each column in
%   DATAMATRIX.
%
%   VIOLINPLOT(TABLE), VIOLINPLOT(STRUCT), VIOLINPLOT(DATASET)
%   plots violins for each column in TABLE, each field in STRUCT, and
%   each variable in DATASET. The violins are labeled according to
%   the table/dataset variable name or the struct field name.
%
%   VIOLINPLOT(DATAMATRIX, CATEGORYNAMES) plots violins for each
%   column in DATAMATRIX and labels them according to the names in the
%   cell-of-strings CATEGORYNAMES.
%
%   VIOLINPLOT(DATA, CATEGORIES) where double vector DATA and vector
%   CATEGORIES are of equal length; plots violins for each category in
%   DATA.
%
%   violins = VIOLINPLOT(...) returns an object array of
%   <a href="matlab:help('Violin')">Violin</a> objects.
%
%   VIOLINPLOT(..., 'PARAM1', val1, 'PARAM2', val2, ...)
%   specifies optional name/value pairs for all violins:
%     'Width'        Width of the violin in axis space.
%                    Defaults to 0.3
%     'Bandwidth'    Bandwidth of the kernel density estimate.
%                    Should be between 10% and 40% of the data range.
%     'ViolinColor'  Fill color of the violin area and data points.
%                    Defaults to the next default color cycle.
%     'ViolinAlpha'  Transparency of the violin area and data points.
%                    Defaults to 0.3.
%     'EdgeColor'    Color of the violin area outline.
%                    Defaults to [0.5 0.5 0.5]
%     'BoxColor'     Color of the box, whiskers, and the outlines of
%                    the median point and the notch indicators.
%                    Defaults to [0.5 0.5 0.5]
%     'MedianColor'  Fill color of the median and notch indicators.
%                    Defaults to [1 1 1]
%     'ShowData'     Whether to show data points.
%                    Defaults to true
%     'ShowNotches'  Whether to show notch indicators.
%                    Defaults to false
%     'ShowMean'     Whether to show mean indicator
%                    Defaults to false

% Copyright (c) 2016, Bastian Bechtold
% This code is released under the terms of the BSD 3-clause license

[num_samples,x_range] = size(data);

for n=1:x_range
    thisData = data(:, n);
    x_pos = x(n);
    violins(n) = Violin(thisData, x_pos,median_flag, varargin{:});
end

