function [wbName,wbH] = myWaitbar(i,N,wbName,flag)
if isempty(wbName)
    wbH = [];
    return
end

switch flag

    case 'init'
        % multiWaitbar('CloseAll');                                           % close previous waitbars if any
        % wbName = getWaitBarTitle(iSet,NSet);
        wbName = 'computing set of solutions...';
        [~,wbH] = multiWaitbar(wbName, 0,...
            'Color', [0.0039,0.3098,0.5255] );
        wbH.Name = 'samwise solver';
        col = [0 0 0]+0.2;
        col2 = [0.6627, 0.8392, 0.8980];
        wbH.Children(1).BackgroundColor = col;
        wbH.Children(2).BackgroundColor = col;
        wbH.Children(3).BackgroundColor = col;
        wbH.Children(4).BackgroundColor = col;
        wbH.Children(5).BackgroundColor = col;
        wbH.Children(3).ForegroundColor = col2;
        wbH.Children(4).ForegroundColor = col2;

    case 'title'

        waitBarNameNew     = getWaitBarTitle(iSet,NSet);
        [~,wbH] = multiWaitbar(wbName,'RELABEL',waitBarNameNew );
        wbName = waitBarNameNew;

    case 'done'
        waitBarNameNew = 'done';
        [~,wbH] = multiWaitbar(wbName,'RELABEL',waitBarNameNew);
        wbName = waitBarNameNew;

    case 'value'

        [~,wbH] = multiWaitbar(wbName,i/N);

    otherwise
        wbH = [];

end



end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
