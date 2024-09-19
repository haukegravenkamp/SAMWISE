function [wbName,wbH] = initializeWaitbar(showWaitbars)

[wbName,wbH] = myWaitbar(0,[],showWaitbars,'init');                         % initialize wait bar if requested
if ~isempty(wbH)
    nWaitBars = round(numel(wbH.Children)/5);
    if nWaitBars > 0
        wbUpdate = ['computing set of solutions ',num2str(nWaitBars)];
        multiWaitbar(wbName,'RELABEL',wbUpdate );
        wbName = wbUpdate;
    end
end
end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
