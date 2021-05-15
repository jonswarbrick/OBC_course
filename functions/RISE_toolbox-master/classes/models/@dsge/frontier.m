% FRONTIER Rolling efficient frontier.
%    [PORTWGTS,ALLMEAN,ALLCOV] = FRONTIER(UNIVERSE,WINDOW,OFFSET,NPORTS,ACTIVEMAP,CONSET,NUMNONNAN)
%    generates a surface of efficient frontiers given data indicating which
%    stocks are active at each date, portfolio constraints and the number of
%    non missing data points needed for each window showing how asset 
%    allocation influences risk and return over time.   
% 
%    [PW,AM,AC] = FRONTIER(UNIVERSE,WINDOW,OFFSET,NPORTS)
%    generates a surface of efficient frontiers showing how asset allocation 
%    influences risk and return over time.   
% 
%    Inputs:
% 
%       UNIVERSE - Portfolio matrix containing the total return data for a
%                  group of securities.  It is an Mx(N+1) matrix where column 1
%                  contains MATLAB date numbers and the remaining columns 
%                  are the total return data for each security.  N is the
%                  number of securities.
%       WINDOW - Number of periods of data to use to calculate each frontier.
%       OFFSET - Increment in number of periods between each frontier.
%       NPORTS - The number of portfolios to calculate on each frontier.
%       
%    Optional inputs:
% 
% 
%       ACTIVEMAP - An MxN matrix with boolean elements that correspond to
%                   the UNIVERSE where each element indicates if the asset
%                   is part of the UNIVERSE on the corresponding date.  The
%                   default map is an MxN matrix of ones indicating that all
%                   assets are active at all dates.   Note if constraints,
%                   CONSET, other than the default values are used,
%                   the default ACTIVEMAP must be used.
%       CONSET - Portfolio constraints.  Default constraints are generated by 
%                PORTCONS('Default',NumAssets).  This single constraint
%                matrix is applied to each frontier.  Note if an ACTIVEMAP
%                other than the default map, the default constraints must be 
%                used. 
%       NUMNONNAN - Minimum number of non NaN points for each active asset 
%                   in each window of data needed to perform the
%                   optimization.  The default value is WINDOW - N.
% 
%    Outputs: 
%   
%      PORTWGTS - A (Number of Curves)x1 cell array where each element is an 
%      NPORTSx(Number of Assets) matrix of weights allocated to each
%      asset.  Number of Assets = length(UNIVERSE). 
% 
%      ALLMEAN - A (Number of Curves)x1 cell array where each element is an 
%      1 x (Number of Assets) vector of the expected asset returns used to
%      generate each curve on the surface.
%             
%      ALLCOV - A (Number of Curves)x1 cell array where each element is an 
%      (Number of Assets) x (Number of Assets) vector of the covariance
%      matrix used to generate each curve on the surface.
% 
%    See also PORTCONS, PORTOPT.
%
%    Reference page in Doc Center
%       doc frontier
%
%    Other functions named frontier
%
%       dsge/frontier
%