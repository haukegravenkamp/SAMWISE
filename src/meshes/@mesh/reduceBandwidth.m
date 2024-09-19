function [K,indb,varargout]=reduceBandwidth(K,varargin)
%% reduce bandwith of the matrix K and sort all other vectors/matrices in varargin accordingly

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
indr = symrcm(K);                                                           % indices to reduce bandwidth
[~,indb]=sort(indr);                                                        % indices to reverse the ordering
K=K(indr,indr);                                                             % K, reduced bandwidth
nout=length(varargin);                                                      % number of input arrays
varargout{nout}=[];
nd=numel(indr);

for i=1:nout

    if (size(varargin{i},1)==nd) && (size(varargin{i},2)==nd)               % if square matrix
        varargout{i}=varargin{i}(indr,indr);                                % sort both dimensions according to indr
    elseif (size(varargin{i},1)==nd)                                        % first direction to be sorted
        varargout{i}=varargin{i}(indr,:);                                   % sort according to indr
    elseif (size(varargin{i},2)==nd)                                        % second direction to be sorted
        varargout{i}=varargin{i}(:,indr);                                   % sort according to indr
    end

end

end
