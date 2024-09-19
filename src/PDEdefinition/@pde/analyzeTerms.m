
% Checks which terms of the weak form are present
function analyzeTerms(obj)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
for iPde=1:numel(obj)
    
    for iTerm = 1:numel(obj(iPde).terms)
        
        existbx1=~isempty(obj(iPde).terms(iTerm).bx1);                      % check which spatial derivatives exist
        existbz1=~isempty(obj(iPde).terms(iTerm).bz1);
        existbx2=~isempty(obj(iPde).terms(iTerm).bx2);
        existbz2=~isempty(obj(iPde).terms(iTerm).bz2);
        
        existFunc=~isempty(obj(iPde).terms(iTerm).func);                    % whether function of spatial coordinates exists
        
        if obj(iPde).terms(iTerm).orderTime==0                              % no time derivative
            j=1;                                                            % stiffness term
        elseif obj(iPde).terms(iTerm).orderTime==1                          % 1st order time derivative
            j=2;                                                            % mass term
        elseif obj(iPde).terms(iTerm).orderTime==2                          % 2nd order time derivative
            j=3;                                                            % damping term
        end
        
        if all(~[existbx1,existbz1,existbx2,existbz2])                      % if term is 0th order in space: v*u
            i=1;                                                            % term 1
        elseif all([~existbx1,~existbz1,existbx2,existbz2])                 % if term is 1st order in space: v*u'
            i=2;                                                            % term 2
        elseif all([existbx1,existbz1,~existbx2,~existbz2])                 % if term is 1st order in space: v'*u
            i=3;                                                            % term 3
        elseif all([existbx1,existbz1,existbx2,existbz2])                   % if term is 2nd order in space
            i=4;                                                            % term 4
        end
        
        obj(iPde).termNumbers(i,j)=iTerm;
        if existFunc
            obj(iPde).termFuncExist(i,j)=existFunc;
            obj(iPde).termFunc{i,j}=obj(iPde).terms(iTerm).func;
        end
    end
end

end
