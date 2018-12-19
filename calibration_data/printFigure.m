function printFigure( filename, sz_paper, varargin )
% PRINTFIGURE  Print a figure to an EPS file at a given size
%
% ## Syntax
% printFigure( filename, sz_paper [, orientation, sz_figure, fg] )
%
% ## Description
% printFigure( filename, sz_paper [, orientation, sz_figure, fg] )
%   Saves a print of the figure to the given filename, with the given paper
%   size and printing options.
%
% ## Input Arguments
%
% filename -- Output file
%   A character vector containing the path of the output file.
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
%   The print size of the figure's axes. A two-element vector containing
%   the width and height, respectively, of the axes in inches. Note that
%   the width and height are in the context of the figure's screen
%   orientation, not its orientation on the page. Defaults to "best fit" if
%   empty, or if not passed.
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
% - The output resolution is 600 dpi.
% - In the past, this function was used to generate PDF files, but the PDF
%   generator failed for large image sizes. While the EPS generator works
%   more reliably, the output files do not have margins once converted to
%   PDF. Therefore, the paper and figure sizes should be set to the same
%   values, as the paper size seems to be irrelevant.
% - Use the `epstopdf` program that comes with a LaTeX distribution to
%   convert the Encapsulated PostScript file to a PDF.
% - Use the `pdfcrop` program that comes with a LaTeX distribution to
%   eliminate the whitespace around a PDF figure, if needed.
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

set(fg, 'PaperUnits', 'inches');
set(fg, 'PaperSize', sz_paper);
set(fg, 'PaperOrientation', orientation);
set(fg, 'PaperPositionMode', 'manual');

if isempty(sz_figure)
    print(fg, filename, '-painters', '-deps', '-r600', '-bestfit');
else
    margin_bound = 0.25;
    if (is_portrait && (sz_paper(1) <= sz_paper(2))) || (~is_portrait && (sz_paper(1) >= sz_paper(2)))
        space = sz_paper - sz_figure;
    else
        space = flip(sz_paper) - sz_figure;
    end
    
    if any(space < (margin_bound * 2))
        warning('Figure size is too large for the page, taking into account a margin of %f inches.', margin_bound)
    end
    set(fg, 'PaperPosition', [margin_bound, margin_bound, sz_figure(1), sz_figure(2)]);
    print(fg, filename, '-painters', '-deps', '-r600');
end

end

