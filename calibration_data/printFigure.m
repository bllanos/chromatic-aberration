function printFigure( filename, sz_paper, varargin )
% PRINTFIGURE  Print a figure to PDF at a given size
%
% ## Syntax
% printFigure( filename, sz_paper [, orientation, sz_figure, fg] )
%
% ## Description
% printFigure( filename, sz_paper [, orientation, sz_figure, fg] )
%   Saves a PDF of the figure to the given filename, with the given paper
%   size and printing options.
%
% ## Input Arguments
%
% filename -- Output file
%   A character vector containing the path of the output PDF file.
%
% sz_paper -- Paper dimensions
%   A two-element vector containing the width and height, respectively, of
%   the paper on which the figure is to be printed, in inches.
%
% orientation -- Figure orientation
%   A character vector, either 'portrait' or 'landscape', specifying the
%   orientation of the figure on the page. Defaults to 'portrait' if empty,
%   or if not passed.
%
% sz_figure -- Figure dimensions
%   The print size of the figure. A two-element vector containing the width
%   and height, respectively, of the figure in inches. An error will be
%   thrown if the print size results in insufficient margins (less than
%   0.25 inches). Note that the width and height are in the context of the
%   figure's screen orientation, not its orientation on the page. Defaults
%   to "best fit" if empty, or if not passed.
%
% fg -- Figure handle
%   The handle to the figure to save. Defaults to the value of `gcf()` if
%   empty, or if not passed.
%
% ## Side-effects
% - Printing-related properties of the figure are modified.
%
% ## Notes
% - The figure will be anchored at the bottom left of the page.
% - The output resolution is 300 dpi.
% - Use the `pdfcrop` program that comes with a LaTeX distribution to
%   eliminate the whitespace around the figure, if needed.
%
% ## References
% - https://www.mathworks.com/matlabcentral/answers/37318-how-to-save-a-figure-that-is-larger-then-the-screen#answer_116509
%
% See also print

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 13, 2018

nargoutchk(0,0)
narginchk(2,5)

orientation = [];
sz_figure = [];
fg = [];
if ~isempty(varargin)
    orientation = varargin{1};
    if length(varargin) > 1
        sz_figure = varargin{2};
        if length(varargin) > 2
            fg = varargin{3};
        end
    end
end
if isempty(orientation)
    orientation = 'portrait';
end
if isempty(fg)
    fg = gcf;
end

if strcmp(orientation, 'portrait')
    is_portrait = true;
elseif strcmp(orientation, 'landscape')
    is_portrait = false;
else
    error('Unknown paper orientation "%s".', orientation);
end

if (is_portrait && (sz_paper(1) <= sz_paper(2))) || (~is_portrait && (sz_paper(1) > sz_paper(2)))
    sz_figure_oriented = sz_figure;
else
    sz_figure_oriented = flip(sz_figure);
end

margin_bound = 0.25;
if any(sz_paper - sz_figure_oriented) < (margin_bound * 2)
    error('Figure size is too large for the page, taking into account a margin of %f inches.', margin_bound)
end

set(fg, 'PaperUnits', 'inches');
set(fg, 'PaperSize', sz_paper);
set(fg, 'PaperOrientation', orientation);
set(fg, 'PaperPositionMode', 'manual');
set(fg, 'PaperPosition', [margin_bound, margin_bound, sz_figure_oriented(1), sz_figure_oriented(2)]);

print(fg, filename, '-dpdf', '-r300');

end

